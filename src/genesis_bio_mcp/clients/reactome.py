"""Reactome pathway analysis client."""

from __future__ import annotations

import asyncio
import logging

import httpx

from genesis_bio_mcp.models import Pathway, PathwayContext

logger = logging.getLogger(__name__)

_REACTOME_ANALYSIS_URL = "https://reactome.org/AnalysisService/identifiers/"
_REACTOME_PATHWAYS_URL = (
    "https://reactome.org/AnalysisService/download/{token}/pathways/TOTAL/result.json"
)

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

    async def get_pathway_context(
        self, gene_symbol: str, max_pathways: int = 10
    ) -> PathwayContext | None:
        """Submit gene to Reactome analysis and return top enriched pathways."""
        async with _SEMAPHORE:
            token = await self._run_analysis(gene_symbol)
            if token is None:
                return None

            pathways = await self._fetch_pathways(token, max_pathways)
            if not pathways:
                return None

            top_name = pathways[0].display_name if pathways else None
            return PathwayContext(
                gene_symbol=gene_symbol,
                pathways=pathways,
                top_pathway_name=top_name,
            )

    async def _run_analysis(self, gene_symbol: str) -> str | None:
        """POST gene identifier to Reactome analysis service, return token."""
        try:
            resp = await self._client.post(
                _REACTOME_ANALYSIS_URL,
                content=gene_symbol,
                params={
                    "interactors": "false",
                    "pageSize": 20,
                    "page": 1,
                    "sortBy": "ENTITIES_PVALUE",
                    "order": "ASC",
                    "resource": "TOTAL",
                    "pValue": 1,
                    "includeDisease": "true",
                    "min": 1,
                    "species": "Homo sapiens",
                },
                headers={"Content-Type": "text/plain"},
                timeout=25.0,
            )
            resp.raise_for_status()
            data = resp.json()
            return data.get("summary", {}).get("token")
        except Exception as exc:
            logger.warning("Reactome analysis failed for %s: %s", gene_symbol, exc)
            return None

    async def _fetch_pathways(self, token: str, max_pathways: int) -> list[Pathway]:
        """Fetch pathway results for a completed Reactome analysis token."""
        try:
            resp = await self._client.get(
                _REACTOME_PATHWAYS_URL.format(token=token),
                timeout=20.0,
            )
            resp.raise_for_status()
            data = resp.json()
            pathways_data = data.get("pathways", []) or []
        except Exception as exc:
            logger.warning("Reactome pathway fetch failed for token %s: %s", token, exc)
            return []

        pathways: list[Pathway] = []
        for p in pathways_data[:max_pathways]:
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
