"""Reactome pathway analysis client."""

from __future__ import annotations

import asyncio
import logging

import httpx

from genesis_bio_mcp.models import Pathway, PathwayContext

logger = logging.getLogger(__name__)

_REACTOME_ANALYSIS_URL = "https://reactome.org/AnalysisService/identifiers/"
_REACTOME_CONTENT_URL = "https://reactome.org/ContentService"

_SEMAPHORE = asyncio.Semaphore(3)

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
        # Session-scoped cache: only successful results are cached.
        # Failures (network errors, empty results) are NOT cached so that transient
        # errors don't permanently poison the cache for the remainder of the session.
        self._cache: dict[str, PathwayContext] = {}

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

        The AnalysisService POST /identifiers/ response includes pathways inline —
        no separate download step is needed. Using the inline result avoids a second
        round-trip and the separate /download/ endpoint (which has had availability
        issues historically).
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
        if not pathways_data:
            token = data.get("summary", {}).get("token")
            logger.debug(
                "Reactome inline pathways empty for %s (token=%s); total found=%s",
                gene_symbol,
                token,
                data.get("pathwaysFound"),
            )
            return []

        return _parse_pathways(pathways_data[:max_pathways])

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

            return await self._fetch_pathway_genes(stid, max_genes)

    async def _search_pathway_stid(self, name: str) -> str | None:
        """Search Reactome ContentService for a pathway by name; return its stId."""
        try:
            resp = await self._client.get(
                f"{_REACTOME_CONTENT_URL}/data/search/query",
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
            resp = await self._client.get(
                f"{_REACTOME_CONTENT_URL}/data/participants/{stid}",
                params={"type": "ReferenceGeneProduct"},
                timeout=20.0,
            )
            resp.raise_for_status()
            data = resp.json()
        except Exception as exc:
            logger.warning("Reactome participant fetch failed for %s: %s", stid, exc)
            return []

        genes: set[str] = set()
        for entry in data:
            # Each entry is a PhysicalEntity; its refEntities carry the gene names
            ref_entities = entry.get("refEntities", [])
            for ref in ref_entities:
                gene_names = ref.get("geneName", [])
                if isinstance(gene_names, list):
                    genes.update(g.upper() for g in gene_names if g)
                elif isinstance(gene_names, str) and gene_names:
                    genes.add(gene_names.upper())

            # Fallback: identifier field on the entry itself
            if not genes:
                identifier = entry.get("identifier", "")
                if identifier:
                    genes.add(identifier.upper())

        return sorted(genes)[:max_genes]


def _parse_pathways(pathways_data: list[dict]) -> list[Pathway]:
    """Parse raw Reactome pathway dicts into Pathway models."""
    pathways: list[Pathway] = []
    for p in pathways_data:
        p_value = p.get("entities", {}).get("pValue")
        gene_count = p.get("entities", {}).get("total")
        name = p.get("name", "")
        reactome_id = p.get("stId", "")
        category = _infer_category(name)

        try:
            p_float = float(p_value) if p_value is not None else None
        except (ValueError, TypeError):
            p_float = None

        pathways.append(
            Pathway(
                reactome_id=reactome_id,
                display_name=name,
                p_value=p_float,
                gene_count=int(gene_count) if gene_count is not None else None,
                category=category,
            )
        )

    return pathways


def _infer_category(pathway_name: str) -> str | None:
    """Infer a broad category from the pathway name using keyword matching."""
    name_lower = pathway_name.lower()
    for keyword, category in _KEYWORD_CATEGORIES:
        if keyword in name_lower:
            return category
    return None
