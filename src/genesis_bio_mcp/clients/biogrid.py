"""BioGRID protein–protein interaction client.

BioGRID (Biological General Repository for Interaction Datasets) provides
curated PPI data from published literature.  Access requires a free API key:
  https://wiki.thebiogrid.org/doku.php/biogrid_rest_api

Set the BIOGRID_ACCESS_KEY environment variable before starting the server.
If the key is absent, all tool calls return a descriptive error message.
"""

from __future__ import annotations

import asyncio
import logging
import os

import httpx

from genesis_bio_mcp.models import BioGRIDInteraction, BioGRIDInteractome

logger = logging.getLogger(__name__)

_BIOGRID_BASE = "https://webservice.thebiogrid.org/interactions/"
_SEMAPHORE = asyncio.Semaphore(2)


class BioGRIDClient:
    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client
        self._cache: dict[str, BioGRIDInteractome] = {}

    async def get_interactions(
        self,
        gene_symbol: str,
        max_results: int = 500,
    ) -> BioGRIDInteractome | None:
        """Return curated PPI interactions for *gene_symbol* from BioGRID.

        Returns ``None`` if the API key is not set.  Returns an empty
        ``BioGRIDInteractome`` (``total_interactions=0``) if BioGRID returns
        no records for the gene.
        """
        access_key = os.environ.get("BIOGRID_ACCESS_KEY")
        if not access_key:
            logger.warning("BIOGRID_ACCESS_KEY not set; BioGRID queries disabled")
            return None

        symbol = gene_symbol.strip().upper()
        if symbol in self._cache:
            logger.debug("BioGRID cache hit: %s", symbol)
            return self._cache[symbol]

        async with _SEMAPHORE:
            result = await self._fetch(symbol, access_key, max_results)

        if result is not None:
            self._cache[symbol] = result
        return result

    async def _fetch(
        self, symbol: str, access_key: str, max_results: int
    ) -> BioGRIDInteractome | None:
        params = {
            "searchNames": "true",
            "geneList": symbol,
            "taxId": "9606",
            "format": "json",
            "accessKey": access_key,
            "includeInteractors": "true",
            "max": str(max_results),
        }
        try:
            resp = await self._client.get(
                _BIOGRID_BASE,
                params=params,
                timeout=25.0,
            )
            if resp.status_code == 403:
                logger.error("BioGRID returned 403 — check BIOGRID_ACCESS_KEY")
                return None
            resp.raise_for_status()
            data: dict = resp.json()
        except Exception as exc:
            logger.warning("BioGRID fetch failed for %s: %s", symbol, exc)
            return None

        interactions: list[BioGRIDInteraction] = []
        partners: set[str] = set()

        for record in data.values():
            if not isinstance(record, dict):
                continue
            gene_a = record.get("OFFICIAL_SYMBOL_A", "")
            gene_b = record.get("OFFICIAL_SYMBOL_B", "")
            partner = gene_b if gene_a.upper() == symbol else gene_a
            partners.add(partner)
            interactions.append(
                BioGRIDInteraction(
                    interactor_a=gene_a,
                    interactor_b=gene_b,
                    experimental_system=record.get("EXPERIMENTAL_SYSTEM"),
                    experimental_system_type=record.get("EXPERIMENTAL_SYSTEM_TYPE"),
                    pubmed_id=str(record["PUBMED_ID"]) if record.get("PUBMED_ID") else None,
                    throughput=record.get("THROUGHPUT"),
                )
            )

        return BioGRIDInteractome(
            gene_symbol=symbol,
            total_interactions=len(interactions),
            unique_partners=len(partners),
            interactions=interactions[:50],
        )
