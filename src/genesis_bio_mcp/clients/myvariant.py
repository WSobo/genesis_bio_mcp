"""MyVariant.info client — consolidated variant-effect lookup.

MyVariant.info (https://myvariant.info, maintained by TSRI/SuLabs) aggregates
ClinVar, dbNSFP (AlphaMissense, REVEL, CADD, SIFT, PolyPhen, ...), and
gnomAD into a single HGVS-keyed REST API. One GET on a genomic HGVS
coordinate returns clinical significance, pathogenicity predictions, and
population frequency without the per-source search fragility of E-utils.

Assembly: we always request ``hg38`` because the rest of the project uses
GRCh38 coordinates (gnomAD v4, Ensembl, etc.).
"""

from __future__ import annotations

import asyncio
import logging
from typing import Any

import httpx

from genesis_bio_mcp.models import (
    ClinVarAssertion,
    ClinVarRecord,
    InSilicoPredictions,
    PopulationFrequency,
    VariantAnnotation,
)

logger = logging.getLogger(__name__)

_BASE_URL = "https://myvariant.info/v1"
_SEMAPHORE = asyncio.Semaphore(3)

# Only pull the fields we actually surface — keeps response size small and
# avoids depending on dbNSFP fields that vary by release.
_FIELDS = (
    "clinvar.hgvs,clinvar.rsid,clinvar.variant_id,clinvar.rcv,"
    "dbnsfp.alphamissense,dbnsfp.revel,dbnsfp.cadd,dbnsfp.sift,dbnsfp.polyphen2,"
    "gnomad_exome.af"
)


class MyVariantClient:
    """Client for the MyVariant.info REST API.

    Caches lookups by HGVS genomic string within the session; a variant's
    annotation is stable for the lifetime of an MCP server run.
    """

    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client
        self._cache: dict[str, VariantAnnotation | None] = {}

    async def get_annotation(self, hgvs_genomic: str) -> VariantAnnotation | None:
        """Fetch consolidated variant annotation by GRCh38 genomic HGVS.

        Args:
            hgvs_genomic: HGVS genomic string, e.g. ``"chr17:g.7675088C>T"``.

        Returns:
            :class:`VariantAnnotation` on success, ``None`` if no record
            exists or the API is unreachable.
        """
        key = hgvs_genomic.strip()
        if key in self._cache:
            logger.debug("MyVariant cache hit: %s", key)
            return self._cache[key]
        async with _SEMAPHORE:
            result = await self._fetch(key)
        self._cache[key] = result
        return result

    async def _fetch(self, hgvs: str) -> VariantAnnotation | None:
        url = f"{_BASE_URL}/variant/{hgvs}"
        try:
            resp = await self._client.get(
                url,
                params={"assembly": "hg38", "fields": _FIELDS},
                timeout=20.0,
            )
            if resp.status_code == 404:
                return None
            resp.raise_for_status()
            data = resp.json()
        except Exception as exc:
            logger.warning("MyVariant fetch failed for %s: %s", hgvs, exc)
            return None
        if not isinstance(data, dict) or not data or data.get("notfound"):
            return None
        return _parse_annotation(hgvs, data)

    async def query_by_protein_change(
        self, gene_symbol: str, position: int, alt_aa: str
    ) -> VariantAnnotation | None:
        """Query MyVariant by gene symbol + protein position + alt amino acid.

        Used when no gnomAD variant_id is available (e.g., somatic cancer
        hotspots absent from the population reference). Returns the first
        matching hit from MyVariant's /query endpoint, or None on miss/error.

        Args:
            gene_symbol: HGNC symbol (e.g. ``"BRAF"``).
            position: 1-indexed protein position (e.g. ``600``).
            alt_aa: Alt amino acid one-letter code (e.g. ``"E"``).
        """
        cache_key = f"_pq:{gene_symbol.upper()}:{position}:{alt_aa.upper()}"
        if cache_key in self._cache:
            logger.debug("MyVariant protein-query cache hit: %s", cache_key)
            return self._cache[cache_key]
        # ClinVar.gene.symbol is the most reliable gene anchor across dbNSFP
        # transcripts; aa.pos + aa.alt restrict to the exact missense change.
        q = (
            f"clinvar.gene.symbol:{gene_symbol.upper()} "
            f"AND dbnsfp.aa.pos:{position} "
            f"AND dbnsfp.aa.alt:{alt_aa.upper()}"
        )
        async with _SEMAPHORE:
            try:
                resp = await self._client.get(
                    f"{_BASE_URL}/query",
                    params={"q": q, "assembly": "hg38", "fields": _FIELDS, "size": 1},
                    timeout=20.0,
                )
                resp.raise_for_status()
                data = resp.json()
            except Exception as exc:
                logger.warning(
                    "MyVariant protein-query failed for %s p.%d?>%s: %s",
                    gene_symbol,
                    position,
                    alt_aa,
                    exc,
                )
                return None
        hits = data.get("hits", []) if isinstance(data, dict) else []
        if not hits:
            return None
        hit = hits[0]
        hgvs = hit.get("_id") or f"{gene_symbol}:p.{position}{alt_aa}"
        result = _parse_annotation(hgvs, hit)
        self._cache[cache_key] = result
        return result


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------


def _parse_annotation(hgvs: str, data: dict[str, Any]) -> VariantAnnotation:
    clinvar = _parse_clinvar(data.get("clinvar"))
    pops = _parse_frequency(data.get("gnomad_exome"))
    in_silico = _parse_in_silico(data.get("dbnsfp"))
    return VariantAnnotation(
        query=hgvs,
        clinvar=clinvar,
        gnomad=pops,
        in_silico=in_silico,
    )


def _parse_clinvar(raw: Any) -> ClinVarRecord | None:
    if not isinstance(raw, dict):
        return None
    rcv_raw = raw.get("rcv")
    if isinstance(rcv_raw, dict):
        rcv_raw = [rcv_raw]
    elif not isinstance(rcv_raw, list):
        rcv_raw = []
    assertions: list[ClinVarAssertion] = []
    for r in rcv_raw:
        if not isinstance(r, dict):
            continue
        conditions = r.get("conditions")
        if isinstance(conditions, dict):
            condition_names = [conditions.get("name")] if conditions.get("name") else []
        elif isinstance(conditions, list):
            condition_names = [
                c.get("name") for c in conditions if isinstance(c, dict) and c.get("name")
            ]
        else:
            condition_names = []
        assertions.append(
            ClinVarAssertion(
                accession=r.get("accession") or "",
                significance=r.get("clinical_significance") or "Unknown",
                review_status=r.get("review_status") or "",
                origin=r.get("origin") or "",
                last_evaluated=r.get("last_evaluated") or None,
                conditions=[c for c in condition_names if c],
            )
        )
    significance = _summarize_significance(assertions)
    return ClinVarRecord(
        rsid=raw.get("rsid") or None,
        hgvs_protein=_first(raw.get("hgvs", {}).get("protein")),
        hgvs_coding=_first(raw.get("hgvs", {}).get("coding")),
        hgvs_genomic=_first(raw.get("hgvs", {}).get("genomic")),
        significance_summary=significance,
        assertions=assertions,
    )


def _first(xs: Any) -> str | None:
    """Return the first string from an iterable of strings, or None."""
    if isinstance(xs, list) and xs:
        s = xs[0]
        return s if isinstance(s, str) else None
    if isinstance(xs, str):
        return xs
    return None


def _summarize_significance(assertions: list[ClinVarAssertion]) -> str:
    """Collapse multiple RCV records into a single significance verdict.

    Precedence (strongest first): Pathogenic → Likely pathogenic →
    Conflicting → Uncertain significance → Likely benign → Benign.
    """
    if not assertions:
        return "No classification"
    order = [
        "Pathogenic",
        "Likely pathogenic",
        "Conflicting classifications of pathogenicity",
        "Conflicting interpretations of pathogenicity",
        "Uncertain significance",
        "Likely benign",
        "Benign",
    ]
    sigs = {a.significance for a in assertions}
    for label in order:
        if label in sigs:
            return label
    return next(iter(sigs))


def _parse_frequency(raw: Any) -> PopulationFrequency | None:
    if not isinstance(raw, dict):
        return None
    af = raw.get("af")
    if isinstance(af, dict):
        overall = af.get("af")
    else:
        overall = af
    if overall is None:
        return None
    try:
        freq = float(overall)
    except (TypeError, ValueError):
        return None
    # Pick a handful of population-stratified fields that gnomAD exposes,
    # safely ignoring missing keys.
    subset: dict[str, float] = {}
    if isinstance(af, dict):
        for key in ("af_afr", "af_amr", "af_eas", "af_nfe", "af_sas"):
            v = af.get(key)
            if isinstance(v, (int, float)):
                subset[key] = float(v)
    return PopulationFrequency(overall_af=freq, by_population=subset)


def _parse_in_silico(raw: Any) -> InSilicoPredictions | None:
    if not isinstance(raw, dict):
        return None
    am = raw.get("alphamissense") or {}
    revel = raw.get("revel") or {}
    cadd = raw.get("cadd") or {}
    sift = raw.get("sift") or {}
    polyphen = raw.get("polyphen2") or {}

    def _score_mean(bundle: dict[str, Any]) -> float | None:
        """Take the mean of per-transcript scores when the field is a list."""
        s = bundle.get("score") if isinstance(bundle, dict) else None
        if isinstance(s, (int, float)):
            return float(s)
        if isinstance(s, list) and s:
            vals = [float(v) for v in s if isinstance(v, (int, float))]
            return sum(vals) / len(vals) if vals else None
        return None

    am_pred = None
    am_pred_raw = am.get("pred") if isinstance(am, dict) else None
    if isinstance(am_pred_raw, list) and am_pred_raw:
        # pred codes: 'P' = likely pathogenic, 'B' = likely benign, 'A' = ambiguous
        majority = max(set(am_pred_raw), key=am_pred_raw.count)
        am_pred = {"P": "likely_pathogenic", "B": "likely_benign", "A": "ambiguous"}.get(majority)
    elif isinstance(am_pred_raw, str):
        am_pred = {"P": "likely_pathogenic", "B": "likely_benign", "A": "ambiguous"}.get(
            am_pred_raw
        )

    return InSilicoPredictions(
        alphamissense_score=_score_mean(am),
        alphamissense_class=am_pred,
        revel_score=_score_mean(revel),
        cadd_phred=cadd.get("phred") if isinstance(cadd.get("phred"), (int, float)) else None,
        sift_score=_score_mean(sift),
        polyphen_score=_score_mean(polyphen),
    )
