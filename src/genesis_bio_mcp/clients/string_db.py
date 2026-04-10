"""STRING protein interaction network client."""

from __future__ import annotations

import asyncio
import logging

import httpx

from genesis_bio_mcp.models import Interactor, ProteinInteractome

logger = logging.getLogger(__name__)

_STRING_NETWORK_URL = "https://string-db.org/api/json/network"
_STRING_RESOLVE_URL = "https://string-db.org/api/json/get_string_ids"

_SEMAPHORE = asyncio.Semaphore(2)

_EVIDENCE_KEYS = [
    "escore",  # experiments
    "dscore",  # databases
    "cscore",  # coexpression
    "tscore",  # textmining
    "hscore",  # homology
    "ascore",  # cooccurrence
    "fscore",  # gene fusion
]

_EVIDENCE_LABELS = {
    "escore": "experiments",
    "dscore": "database",
    "cscore": "coexpression",
    "tscore": "textmining",
    "hscore": "homology",
    "ascore": "cooccurrence",
    "fscore": "gene_fusion",
}


class StringDbClient:
    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client

    async def get_interactome(
        self, gene_symbol: str, required_score: int = 700, limit: int = 20
    ) -> ProteinInteractome | None:
        """Return top interaction partners for a gene from STRING."""
        async with _SEMAPHORE:
            string_id = await self._resolve_string_id(gene_symbol)
            if string_id is None:
                logger.info("STRING: could not resolve '%s' to a STRING ID", gene_symbol)
                return None

            interactions = await self._fetch_network(string_id, required_score, limit)
            if not interactions:
                return ProteinInteractome(
                    gene_symbol=gene_symbol,
                    total_partners=0,
                    top_interactors=[],
                )

            interactors = _parse_interactions(gene_symbol, string_id, interactions)
            return ProteinInteractome(
                gene_symbol=gene_symbol,
                total_partners=len(interactors),
                top_interactors=interactors,
            )

    async def _resolve_string_id(self, gene_symbol: str) -> str | None:
        """Resolve HGNC symbol to STRING ID (9606.ENSP...)."""
        try:
            resp = await self._client.get(
                _STRING_RESOLVE_URL,
                params={
                    "identifiers": gene_symbol,
                    "species": 9606,
                    "limit": 1,
                    "echo_query": 1,
                    "caller_identity": "genesis-bio-mcp",
                },
                timeout=15.0,
            )
            resp.raise_for_status()
            data = resp.json()
            if data:
                return data[0].get("stringId")
            return None
        except Exception as exc:
            logger.warning("STRING resolve failed for %s: %s", gene_symbol, exc)
            return None

    async def _fetch_network(self, string_id: str, required_score: int, limit: int) -> list[dict]:
        """Fetch interaction network edges for a STRING ID."""
        try:
            resp = await self._client.get(
                _STRING_NETWORK_URL,
                params={
                    "identifiers": string_id,
                    "species": 9606,
                    "required_score": required_score,
                    "limit": limit,
                    "caller_identity": "genesis-bio-mcp",
                },
                timeout=20.0,
            )
            resp.raise_for_status()
            return resp.json() or []
        except Exception as exc:
            logger.warning("STRING network fetch failed for %s: %s", string_id, exc)
            return []


def _parse_interactions(
    gene_symbol: str, string_id: str, interactions: list[dict]
) -> list[Interactor]:
    """Parse STRING network edges into Interactor list, excluding self-edges."""
    interactors: dict[str, Interactor] = {}
    symbol_upper = gene_symbol.upper()

    for edge in interactions:
        # Each edge has stringId_A, preferredName_A, stringId_B, preferredName_B, score
        for side in ("A", "B"):
            sid = edge.get(f"stringId_{side}", "")
            name = edge.get(f"preferredName_{side}", "")
            if sid == string_id or name.upper() == symbol_upper:
                continue  # skip the query gene itself

            score = float(edge.get("score", 0)) / 1000.0  # STRING returns 0-1000

            evidence: list[str] = []
            for key, label in _EVIDENCE_LABELS.items():
                val = edge.get(key, 0)
                if isinstance(val, (int, float)) and float(val) > 150:
                    evidence.append(label)

            if name not in interactors or interactors[name].score < score:
                interactors[name] = Interactor(
                    gene_symbol=name,
                    protein_name=name,
                    score=round(score, 3),
                    evidence_types=evidence,
                )

    return sorted(interactors.values(), key=lambda x: x.score, reverse=True)
