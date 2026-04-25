"""Reactome pathway analysis client."""

from __future__ import annotations

import asyncio
import logging

import httpx

from genesis_bio_mcp.config.settings import settings
from genesis_bio_mcp.models import Pathway, PathwayContext

logger = logging.getLogger(__name__)

_REACTOME_ANALYSIS_URL = "https://reactome.org/AnalysisService/identifiers/"
_REACTOME_CONTENT_URL = "https://reactome.org/ContentService"

_SEMAPHORE = asyncio.Semaphore(settings.reactome_semaphore_limit)

# Top-level Reactome categories for human pathways (R-HSA prefixes)
_CATEGORY_MAP = {
    "R-HSA-162582": "Signal Transduction",
    "R-HSA-1430728": "Metabolism",
    "R-HSA-168256": "Immune System",
    "R-HSA-1640170": "Cell Cycle",
    "R-HSA-5633007": "Programmed Cell Death",
    "R-HSA-397014": "Muscle Contraction",
    "R-HSA-212436": "Generic Transcription",
    "R-HSA-9612973": "Autophagy",
    "R-HSA-382551": "Transport of Small Molecules",
    "R-HSA-1474244": "Extracellular Matrix",
}

# Keyword-based category inference for pathway names
_KEYWORD_CATEGORIES = [
    ("signal", "Signaling"),
    ("mapk", "Signaling"),
    ("ras", "Signaling"),
    ("pi3k", "Signaling"),
    ("jak", "Signaling"),
    ("wnt", "Signaling"),
    ("notch", "Signaling"),
    ("hedgehog", "Signaling"),
    ("metabolism", "Metabolism"),
    ("metabolic", "Metabolism"),
    ("biosynthesis", "Metabolism"),
    ("immune", "Immune"),
    ("innate", "Immune"),
    ("adaptive", "Immune"),
    ("interferon", "Immune"),
    ("interleukin", "Immune"),
    ("cell cycle", "Cell Cycle"),
    ("dna repair", "DNA Repair"),
    ("dna replication", "DNA Repair"),
    ("apoptosis", "Apoptosis"),
    ("autophagy", "Autophagy"),
    ("transcription", "Transcription"),
    ("translation", "Translation"),
    ("rna", "RNA Processing"),
    ("splicing", "RNA Processing"),
    ("transport", "Transport"),
    ("vesicle", "Transport"),
    ("extracellular", "ECM"),
    ("collagen", "ECM"),
    ("development", "Development"),
    ("differentiation", "Development"),
]


class ReactomeClient:
    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client
        # Session-scoped caches: only successful results are cached.
        # Failures (network errors, empty results) are NOT cached so that transient
        # errors don't permanently poison the cache for the remainder of the session.
        self._cache: dict[str, PathwayContext] = {}
        self._members_cache: dict[str, list[str]] = {}

    async def get_pathway_context(
        self, gene_symbol: str, max_pathways: int = 10
    ) -> PathwayContext | None:
        """Submit gene to Reactome analysis and return top enriched pathways."""
        symbol = gene_symbol.strip().upper()
        if symbol in self._cache:
            logger.debug("Reactome cache hit: %s", symbol)
            return self._cache[symbol]

        async with _SEMAPHORE:
            pathways = await self._run_analysis(symbol, max_pathways)
            if not pathways:
                # Do NOT cache failures — transient errors must be retryable
                return None

            top_name = pathways[0].display_name if pathways else None
            result = PathwayContext(
                gene_symbol=symbol,
                pathways=pathways,
                top_pathway_name=top_name,
            )
            self._cache[symbol] = result
            return result

    async def _run_analysis(self, gene_symbol: str, max_pathways: int) -> list[Pathway]:
        """POST gene identifier to Reactome AnalysisService; parse inline pathways.

        Attempts to read pathways inline from the POST response first. If the inline
        ``pathways`` array is empty (Reactome sometimes omits it even when
        ``pathwaysFound`` > 0), falls back to fetching via the token endpoint:
        ``GET /AnalysisService/token/{token}/pathways/TOTAL/``.
        """
        try:
            resp = await self._client.post(
                _REACTOME_ANALYSIS_URL,
                # Reactome expects identifiers as newline-separated plain text.
                # A trailing newline ensures the identifier is parsed correctly.
                content=f"{gene_symbol}\n",
                params={
                    "interactors": "false",
                    "pageSize": max_pathways,
                    "page": 1,
                    "sortBy": "ENTITIES_PVALUE",
                    "order": "ASC",
                    "resource": "TOTAL",
                    "pValue": 1,
                    "includeDisease": "true",
                    "min": 1,
                },
                headers={"Content-Type": "text/plain"},
                timeout=30.0,
            )
            resp.raise_for_status()
            data = resp.json()
        except Exception as exc:
            logger.warning("Reactome analysis failed for %s: %s", gene_symbol, exc)
            return []

        pathways_data: list[dict] = data.get("pathways", []) or []
        if pathways_data:
            return _parse_pathways(pathways_data[:max_pathways])

        # Inline pathways absent — fall back to token-based fetch
        token = data.get("summary", {}).get("token")
        logger.warning(
            "Reactome inline pathways empty for %s (token=%s, pathwaysFound=%s); "
            "falling back to token endpoint",
            gene_symbol,
            token,
            data.get("pathwaysFound"),
        )
        if not token:
            return []

        return await self._fetch_pathways_by_token(token, max_pathways)

    async def _fetch_pathways_by_token(self, token: str, max_pathways: int) -> list[Pathway]:
        """GET pathways for a completed Reactome analysis token."""
        try:
            resp = await self._client.get(
                f"https://reactome.org/AnalysisService/token/{token}/pathways/TOTAL/",
                params={
                    "pageSize": max_pathways,
                    "page": 1,
                    "sortBy": "ENTITIES_PVALUE",
                    "order": "ASC",
                    "resource": "TOTAL",
                },
                timeout=20.0,
            )
            resp.raise_for_status()
            # Response is a JSON array of pathway objects
            raw: list[dict] = resp.json()
        except Exception as exc:
            logger.warning("Reactome token fetch failed (token=%s): %s", token, exc)
            return []

        if not isinstance(raw, list):
            # Some Reactome responses wrap the array; handle both forms
            raw = raw.get("pathways", []) if isinstance(raw, dict) else []

        return _parse_pathways(raw[:max_pathways])

    async def get_pathway_members(self, pathway_name_or_id: str, max_genes: int = 50) -> list[str]:
        """Return HGNC gene symbols for all members of a named Reactome pathway.

        Searches ContentService for the pathway by name or stable ID, then fetches
        all ReferenceGeneProduct participants.

        Args:
            pathway_name_or_id: Pathway display name (e.g. "MAPK signaling") or
                Reactome stable ID (e.g. "R-HSA-5673001").
            max_genes: Maximum number of gene symbols to return (default 50).

        Returns:
            Sorted list of unique HGNC gene symbols, or an empty list if the
            pathway cannot be found or has no gene participants.
        """
        async with _SEMAPHORE:
            stid = pathway_name_or_id.strip()

            # If it looks like a stable ID (R-HSA-...) skip the search step
            if not stid.startswith("R-"):
                stid = await self._search_pathway_stid(stid)
                if stid is None:
                    return []

            cache_key = f"{stid}:{max_genes}"
            cached = self._members_cache.get(cache_key)
            if cached is not None:
                logger.debug("Reactome members cache hit: %s", cache_key)
                return cached

            genes = await self._fetch_pathway_genes(stid, max_genes)
            if genes:
                # Only cache successful fetches; transient failures must be retryable.
                self._members_cache[cache_key] = genes
            return genes

    async def _search_pathway_stid(self, name: str) -> str | None:
        """Search Reactome ContentService for a pathway by name; return its stId."""
        try:
            # NOTE: the search endpoint lives at /search/query (no /data/ prefix);
            # /data/search/query 404s.
            resp = await self._client.get(
                f"{_REACTOME_CONTENT_URL}/search/query",
                params={
                    "query": name,
                    "types": "Pathway",
                    "species": "Homo sapiens",
                    "cluster": "true",
                    "rows": 5,
                },
                timeout=15.0,
            )
            resp.raise_for_status()
            data = resp.json()
        except Exception as exc:
            logger.warning("Reactome pathway search failed for '%s': %s", name, exc)
            return None

        results = data.get("results", [])
        for group in results:
            entries = group.get("entries", [])
            if entries:
                return entries[0].get("stId")
        return None

    async def _fetch_pathway_genes(self, stid: str, max_genes: int) -> list[str]:
        """Fetch ReferenceGeneProduct participants for a Reactome pathway stId."""
        try:
            # 60s timeout: large pathways like MAPK signaling return ~300 KB of
            # nested participant JSON; 20s left no margin for parse latency.
            resp = await self._client.get(
                f"{_REACTOME_CONTENT_URL}/data/participants/{stid}",
                params={"type": "ReferenceGeneProduct"},
                timeout=60.0,
            )
            resp.raise_for_status()
            data = resp.json()
        except Exception as exc:
            logger.warning("Reactome participant fetch failed for %s: %s", stid, exc)
            return []

        genes: set[str] = set()
        for entry in data:
            for ref in entry.get("refEntities", []):
                # Only ReferenceGeneProduct entries carry gene-symbol info; the
                # ?type= filter applies to the top-level participant query, not
                # to refEntities, so non-protein entries (e.g. ReferenceMolecule
                # for water/ATP) leak through.
                if ref.get("schemaClass") != "ReferenceGeneProduct":
                    continue
                # displayName format: "UniProt:P12345 SYMBOL"; the gene symbol
                # is the trailing token. There is no "geneName" field on these
                # responses (despite older docs).
                display = ref.get("displayName", "")
                if " " in display:
                    sym = display.rsplit(" ", 1)[-1].strip().upper()
                    if sym:
                        genes.add(sym)
            if len(genes) >= max_genes:
                break

        return sorted(genes)[:max_genes]


def _parse_pathways(pathways_data: list[dict]) -> list[Pathway]:
    """Parse raw Reactome pathway dicts into Pathway models.

    Two-pass dedup:
    1. By Reactome stable ID — the analysis API can return the same stId twice
       across different resource subsets.
    2. By display_name (case-insensitive) — Reactome sometimes returns a parent
       and a child pathway with identical human-readable names but different
       stable IDs. Keeps the row with the smallest p_value.

    Upstream sortBy=ENTITIES_PVALUE puts the strongest hit first, so keeping
    the first occurrence in the stId pass preserves the most informative row.
    """
    seen_ids: set[str] = set()
    parsed: list[Pathway] = []
    for p in pathways_data:
        reactome_id = p.get("stId", "")
        # Bug R (v0.3.4): the AnalysisService doesn't have a species filter
        # at the request level, so it can return pathways from any species
        # the gene maps to (R-CFA-* canine, R-MMU-* mouse, etc.). For human
        # gene queries we only want human pathways. Drop anything that
        # isn't an R-HSA-* stable ID — a pathway with no stId at all is
        # also dropped because it's not actionable.
        if not reactome_id or not reactome_id.startswith("R-HSA-"):
            continue
        if reactome_id in seen_ids:
            continue
        seen_ids.add(reactome_id)

        p_value = p.get("entities", {}).get("pValue")
        gene_count = p.get("entities", {}).get("total")
        name = p.get("name", "")
        category = _infer_category(name)

        try:
            p_float = float(p_value) if p_value is not None else None
        except (ValueError, TypeError):
            p_float = None

        parsed.append(
            Pathway(
                reactome_id=reactome_id,
                display_name=name,
                p_value=p_float,
                gene_count=int(gene_count) if gene_count is not None else None,
                category=category,
            )
        )

    # Name-based dedup: collapse parent/child pathways that share a display name,
    # keeping the row with the smallest p-value.
    by_name: dict[str, Pathway] = {}
    ordering: list[str] = []
    for p in parsed:
        key = p.display_name.strip().lower()
        if not key:
            sentinel = f"__anon_{len(ordering)}__"
            by_name[sentinel] = p
            ordering.append(sentinel)
            continue
        existing = by_name.get(key)
        if existing is None:
            by_name[key] = p
            ordering.append(key)
            continue
        ep = existing.p_value if existing.p_value is not None else float("inf")
        np_ = p.p_value if p.p_value is not None else float("inf")
        if np_ < ep:
            by_name[key] = p

    return [by_name[k] for k in ordering]


def _infer_category(pathway_name: str) -> str | None:
    """Infer a broad category from the pathway name using keyword matching."""
    name_lower = pathway_name.lower()
    for keyword, category in _KEYWORD_CATEGORIES:
        if keyword in name_lower:
            return category
    return None
