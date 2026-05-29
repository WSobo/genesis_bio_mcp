"""Ensembl REST client — gene coordinates + VEP consequence prediction.

Wraps https://rest.ensembl.org. No API key required. Rate limit is 15 req/s
burst and 55 000/hr sustained; a semaphore of 5 keeps us well inside both.

Endpoints used:
- ``GET /lookup/symbol/homo_sapiens/{symbol}`` — gene → Ensembl ID + coords
- ``GET /overlap/region/human/{region}`` — regulatory feature overlap
- ``GET /vep/human/hgvs/{hgvs}`` — consequence prediction from HGVS
- ``GET /vep/human/region/{region}/{allele}`` — consequence by coordinate

VEP responses include every transcript × every consequence, so by default
we filter to the canonical transcript and surface the most-severe term.
Callers can pass ``include_all_transcripts=True`` to get the full list.
"""

from __future__ import annotations

import asyncio
import logging

import httpx

from genesis_bio_mcp.models import (
    EnsemblGene,
    TranscriptInfo,
    VEPConsequence,
    VEPConsequenceReport,
)
from genesis_bio_mcp.tools.variant_parser import canonical_three_letter, parse_protein_change

logger = logging.getLogger(__name__)

_ENSEMBL_BASE = "https://rest.ensembl.org"
_HEADERS = {"Accept": "application/json"}
_SEMAPHORE = asyncio.Semaphore(5)

# Transient Ensembl VEP 5xx errors are common; retry with exponential
# backoff so a one-off flake doesn't poison the aggregated variant report.
_MAX_RETRIES = 3
_BASE_BACKOFF_SECS = 0.5


class EnsemblClient:
    """Session-scoped Ensembl REST client.

    - Gene lookups cached by HGNC symbol (uppercased)
    - VEP results cached by HGVS string (or coordinate tuple)
    - Returns ``None`` on any failure; never raises to caller
    """

    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client
        self._gene_cache: dict[str, EnsemblGene] = {}
        self._vep_cache: dict[str, VEPConsequenceReport] = {}

    async def lookup_gene(self, gene_symbol: str) -> EnsemblGene | None:
        """Resolve an HGNC symbol to an Ensembl gene record with coordinates."""
        symbol = gene_symbol.strip().upper()
        if symbol in self._gene_cache:
            logger.debug("Ensembl gene cache hit: %s", symbol)
            return self._gene_cache[symbol]

        url = f"{_ENSEMBL_BASE}/lookup/symbol/homo_sapiens/{symbol}"
        params = {"expand": "1"}
        async with _SEMAPHORE:
            data = await self._get_json(url, params=params)
        if not data or not isinstance(data, dict):
            return None

        gene = _parse_ensembl_gene(data, symbol)
        if gene is not None:
            self._gene_cache[symbol] = gene
        return gene

    async def get_vep_by_hgvs(
        self,
        hgvs: str,
        *,
        include_all_transcripts: bool = False,
    ) -> VEPConsequenceReport | None:
        """Predict variant consequences for an HGVS string.

        Accepts HGVS.p (with transcript prefix), HGVS.c, or HGVS.g. The
        canonical transcript result is selected by default.
        """
        key = f"hgvs::{hgvs}::all={include_all_transcripts}"
        if key in self._vep_cache:
            return self._vep_cache[key]

        url = f"{_ENSEMBL_BASE}/vep/human/hgvs/{hgvs}"
        async with _SEMAPHORE:
            data = await self._get_json(url)
        if not data or not isinstance(data, list) or not data:
            return None

        report = _parse_vep_response(
            data[0], input_label=hgvs, include_all_transcripts=include_all_transcripts
        )
        if report is not None:
            self._vep_cache[key] = report
        return report

    async def get_vep_by_region(
        self,
        region: str,
        allele: str,
        *,
        include_all_transcripts: bool = False,
    ) -> VEPConsequenceReport | None:
        """Predict consequences from a genomic coordinate region + alt allele.

        ``region`` is the Ensembl region string, e.g. ``"7:140753336-140753336:1"``.
        ``allele`` is the alternate allele (single base for SNVs).
        """
        key = f"region::{region}::{allele}::all={include_all_transcripts}"
        if key in self._vep_cache:
            return self._vep_cache[key]

        url = f"{_ENSEMBL_BASE}/vep/human/region/{region}/{allele}"
        async with _SEMAPHORE:
            data = await self._get_json(url)
        if not data or not isinstance(data, list) or not data:
            return None

        report = _parse_vep_response(
            data[0],
            input_label=f"{region} {allele}",
            include_all_transcripts=include_all_transcripts,
        )
        if report is not None:
            self._vep_cache[key] = report
        return report

    async def get_vep_consequences(
        self,
        gene_symbol: str,
        mutation: str,
        *,
        include_all_transcripts: bool = False,
    ) -> VEPConsequenceReport | None:
        """Convenience: resolve gene → canonical transcript → VEP by HGVS.p.

        Ensembl's ``/vep/human/hgvs`` endpoint rejects bare ``p.Arg175His``
        strings; it requires a transcript accession prefix
        (``ENST00000###.#:p.Arg175His``). We resolve the canonical
        transcript through the gene lookup step and build the query from it.
        """
        gene = await self.lookup_gene(gene_symbol)
        if gene is None or not gene.canonical_transcript_id:
            return None
        # VEP's HGVS endpoint requires canonical protein HGVS
        # (``ENST…:p.Val600Glu``); a bare ``V600E`` is rejected and yields no
        # consequences. Normalize any accepted shorthand ("V600E", "p.V600E",
        # "Val600Glu") — idempotent for callers that already pass "p.Val600Glu".
        try:
            hgvs_p = canonical_three_letter(*parse_protein_change(mutation))
        except ValueError:
            logger.warning(
                "VEP: unparsable protein change '%s' for %s — cannot build HGVS",
                mutation,
                gene_symbol,
            )
            return None
        hgvs = f"{gene.canonical_transcript_id}:{hgvs_p}"
        return await self.get_vep_by_hgvs(hgvs, include_all_transcripts=include_all_transcripts)

    async def _get_json(self, url: str, params: dict | None = None) -> list | dict | None:
        """Issue a GET with retry on 5xx; return parsed JSON or ``None``."""
        backoff = _BASE_BACKOFF_SECS
        for attempt in range(_MAX_RETRIES):
            try:
                resp = await self._client.get(url, params=params, headers=_HEADERS, timeout=25.0)
                if resp.status_code >= 500 and attempt < _MAX_RETRIES - 1:
                    logger.debug(
                        "Ensembl %s returned %s (attempt %s) — retrying",
                        url,
                        resp.status_code,
                        attempt + 1,
                    )
                    await asyncio.sleep(backoff)
                    backoff *= 2
                    continue
                if resp.status_code == 404:
                    return None
                resp.raise_for_status()
                return resp.json()
            except Exception as exc:
                if attempt < _MAX_RETRIES - 1:
                    logger.debug("Ensembl %s failed (%s) — retrying", url, exc)
                    await asyncio.sleep(backoff)
                    backoff *= 2
                    continue
                logger.warning("Ensembl request failed for %s: %s", url, exc)
                return None
        return None


def _parse_ensembl_gene(data: dict, fallback_symbol: str) -> EnsemblGene | None:
    ensembl_id = data.get("id")
    if not ensembl_id:
        return None
    symbol = data.get("display_name") or fallback_symbol
    transcripts_raw = data.get("Transcript", []) or []
    transcripts: list[TranscriptInfo] = []
    canonical_id: str | None = None
    for t in transcripts_raw:
        tid = t.get("id")
        if not tid:
            continue
        is_canonical = bool(t.get("is_canonical"))
        if is_canonical:
            canonical_id = tid
        transcripts.append(
            TranscriptInfo(
                transcript_id=tid,
                is_canonical=is_canonical,
                biotype=t.get("biotype"),
                length=t.get("length"),
            )
        )
    return EnsemblGene(
        ensembl_id=ensembl_id,
        symbol=symbol.upper(),
        chrom=str(data.get("seq_region_name", "")),
        start=int(data["start"]) if data.get("start") is not None else 0,
        end=int(data["end"]) if data.get("end") is not None else 0,
        strand=int(data.get("strand") or 0),
        biotype=data.get("biotype"),
        canonical_transcript_id=canonical_id,
        transcripts=transcripts,
    )


def _parse_vep_response(
    data: dict, *, input_label: str, include_all_transcripts: bool
) -> VEPConsequenceReport | None:
    transcript_consequences = data.get("transcript_consequences") or []
    regulatory = data.get("regulatory_feature_consequences") or []

    parsed: list[VEPConsequence] = []
    for tc in transcript_consequences:
        is_canonical = bool(tc.get("canonical"))
        if not include_all_transcripts and not is_canonical:
            continue
        terms = tc.get("consequence_terms") or []
        parsed.append(
            VEPConsequence(
                consequence_term=", ".join(terms) if terms else "unknown",
                impact=tc.get("impact"),
                transcript_id=tc.get("transcript_id"),
                gene_symbol=tc.get("gene_symbol"),
                biotype=tc.get("biotype"),
                canonical=is_canonical,
                sift_score=tc.get("sift_score"),
                sift_prediction=tc.get("sift_prediction"),
                polyphen_score=tc.get("polyphen_score"),
                polyphen_prediction=tc.get("polyphen_prediction"),
                amino_acids=tc.get("amino_acids"),
                codons=tc.get("codons"),
            )
        )

    regulatory_overlaps: list[str] = []
    for reg in regulatory:
        feature_type = reg.get("biotype") or reg.get("feature_type")
        if feature_type:
            regulatory_overlaps.append(str(feature_type))

    return VEPConsequenceReport(
        input_label=input_label,
        most_severe_consequence=data.get("most_severe_consequence"),
        assembly_name=data.get("assembly_name"),
        consequences=parsed,
        regulatory_overlaps=regulatory_overlaps,
    )
