"""MaveDB deep mutational scanning (DMS) client.

MaveDB (https://www.mavedb.org) is the canonical repository for DMS experiments.
Each "score set" maps thousands of single-amino-acid substitutions or indels to a
quantitative fitness/function score measured in a defined cellular assay.  When a
DMS dataset exists for a target protein it provides the highest-resolution
residue-level tolerance-to-mutation signal available — superior to evolutionary
constraint metrics for predicting the functional consequence of any single variant.

This client searches MaveDB by gene symbol using the text search endpoint and
returns available score-set metadata (URNs, variant counts, publication references).
Individual variant-level scores can be retrieved from MaveDB using the returned URNs.

No API key required.
"""

from __future__ import annotations

import asyncio
import csv
import io
import logging

import httpx

from genesis_bio_mcp.models import DMSResults, DMSScoreSet, MaveDBVariantScore

logger = logging.getLogger(__name__)

_MAVEDB_BASE = "https://api.mavedb.org/api/v1"
_SEARCH_URL = f"{_MAVEDB_BASE}/score-sets/search"
_SEMAPHORE = asyncio.Semaphore(3)


class MaveDBClient:
    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client
        self._cache: dict[str, DMSResults] = {}
        # Per-score-set scores cache: the CSV can be hundreds of kB, so we
        # cache it once per URN rather than re-downloading for every variant.
        self._scores_cache: dict[str, list[dict[str, str]]] = {}

    async def get_dms_scores(self, gene_symbol: str) -> DMSResults | None:
        """Search MaveDB for DMS score sets associated with a gene.

        Uses text search (POST /score-sets/search) which matches against score set
        titles and descriptions.  Returns metadata for all matching score sets sorted
        by variant count descending.  Returns an empty DMSResults (not None) when no
        datasets are found — the absence of DMS data is itself informative.

        Args:
            gene_symbol: HGNC gene symbol, e.g. ``'BRCA1'``.

        Returns:
            :class:`DMSResults` or ``None`` on network error.
        """
        symbol = gene_symbol.strip().upper()
        if symbol in self._cache:
            logger.debug("MaveDB cache hit: %s", symbol)
            return self._cache[symbol]

        async with _SEMAPHORE:
            result = await self._fetch(symbol)

        if result is not None:
            self._cache[symbol] = result
        return result

    async def _fetch(self, symbol: str) -> DMSResults | None:
        payload = {"text": symbol}
        try:
            resp = await self._client.post(_SEARCH_URL, json=payload, timeout=25.0)
            resp.raise_for_status()
            data = resp.json()
        except Exception as exc:
            logger.warning("MaveDB fetch failed for '%s': %s", symbol, exc)
            return None

        if not isinstance(data, dict):
            logger.warning("MaveDB unexpected response type for '%s': %s", symbol, type(data))
            return None

        items: list[dict] = data.get("scoreSets", [])

        if not items:
            return DMSResults(
                gene_symbol=symbol,
                total_score_sets=0,
                total_variants=0,
                score_sets=[],
            )

        score_sets: list[DMSScoreSet] = []
        for item in items:
            urn = item.get("urn") or ""
            if not urn:
                continue

            title = item.get("title") or urn
            short_desc = item.get("shortDescription") or None
            num_variants = int(item.get("numVariants") or 0)
            published_date = item.get("publishedDate") or None

            # Extract target gene / UniProt from targetGenes list
            target_genes: list[dict] = item.get("targetGenes") or []
            target_gene_sym: str | None = None
            uniprot_acc: str | None = None
            if target_genes:
                first = target_genes[0]
                target_gene_sym = first.get("name") or None
                uniprot_acc = first.get("uniprotIdFromMappedMetadata") or None

            # Extract PMID / DOI from primaryPublicationIdentifiers
            pmid: str | None = None
            doi: str | None = None
            for pub in item.get("primaryPublicationIdentifiers") or []:
                db_name = (pub.get("dbName") or "").lower()
                identifier = pub.get("identifier") or ""
                if db_name == "pubmed" and not pmid:
                    pmid = str(identifier)
                elif db_name == "doi" and not doi:
                    doi = str(identifier)

            score_sets.append(
                DMSScoreSet(
                    urn=urn,
                    title=title,
                    short_description=short_desc,
                    num_variants=num_variants,
                    target_gene=target_gene_sym,
                    uniprot_accession=uniprot_acc,
                    published_date=published_date,
                    pmid=pmid,
                    doi=doi,
                )
            )

        # Sort by variant count descending — larger sets = more complete coverage
        score_sets.sort(key=lambda s: s.num_variants, reverse=True)
        total_variants = sum(s.num_variants for s in score_sets)

        return DMSResults(
            gene_symbol=symbol,
            total_score_sets=len(score_sets),
            total_variants=total_variants,
            score_sets=score_sets,
        )

    async def get_variant_score(
        self,
        urn: str,
        hgvs_pro: str,
        score_set_title: str | None = None,
    ) -> list[MaveDBVariantScore]:
        """Return per-variant DMS scores for one score-set matching *hgvs_pro*.

        MaveDB's ``/score-sets/{urn}/scores`` endpoint returns the full CSV of
        every replicate / condition row in the score set. We cache the CSV
        on first access and filter in-memory.

        Args:
            urn: Score-set URN, e.g. ``"urn:mavedb:00000001-a-1"``.
            hgvs_pro: HGVS protein notation to match, e.g. ``"p.Arg175His"``.
            score_set_title: Optional title to attach to each returned score
                (normally populated by the aggregator so the output table
                doesn't force a separate lookup).

        Returns:
            List of :class:`MaveDBVariantScore`, one per matching row. Empty
            list if the variant is not in this score set or the fetch fails.
        """
        rows = await self._load_scores(urn)
        if not rows:
            return []
        matches: list[MaveDBVariantScore] = []
        target = hgvs_pro.strip()
        for row in rows:
            if row.get("hgvs_pro", "").strip() != target:
                continue
            raw_score = row.get("score", "").strip()
            if not raw_score or raw_score == "NA":
                continue
            try:
                score = float(raw_score)
            except ValueError:
                continue
            matches.append(
                MaveDBVariantScore(
                    urn=urn,
                    title=score_set_title or urn,
                    hgvs_pro=target,
                    score=score,
                    epsilon=None,
                )
            )
        return matches

    async def _load_scores(self, urn: str) -> list[dict[str, str]]:
        if urn in self._scores_cache:
            logger.debug("MaveDB scores cache hit: %s", urn)
            return self._scores_cache[urn]
        async with _SEMAPHORE:
            url = f"{_MAVEDB_BASE}/score-sets/{urn}/scores"
            try:
                resp = await self._client.get(url, timeout=30.0)
                resp.raise_for_status()
                text = resp.text
            except Exception as exc:
                logger.warning("MaveDB scores fetch failed for %s: %s", urn, exc)
                self._scores_cache[urn] = []
                return []
        try:
            rows = [dict(row) for row in csv.DictReader(io.StringIO(text))]
        except Exception as exc:
            logger.warning("MaveDB scores CSV parse failed for %s: %s", urn, exc)
            rows = []
        self._scores_cache[urn] = rows
        return rows
