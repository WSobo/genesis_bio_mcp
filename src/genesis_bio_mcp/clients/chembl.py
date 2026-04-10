"""ChEMBL bioactivity client.

Provides quantitative IC50/Ki/Kd potency data from the EMBL-EBI ChEMBL database.
No API key required; rate limit is ~1 req/sec (conservative semaphore at 2).

Flow per gene:
  1. Target search by gene symbol → extract first human single-protein target_chembl_id
  2. Activity query filtered to IC50/Ki/Kd with a pChEMBL value → top 20 by potency
"""

from __future__ import annotations

import asyncio
import logging
from typing import Optional

import httpx

from genesis_bio_mcp.models import ChEMBLActivity, ChEMBLCompounds

logger = logging.getLogger(__name__)

_BASE = "https://www.ebi.ac.uk/chembl/api/data"
_TARGET_SEARCH_URL = f"{_BASE}/target/search"
_ACTIVITY_URL = f"{_BASE}/activity"

_SEMAPHORE = asyncio.Semaphore(2)


class ChEMBLClient:
    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client

    async def get_compounds(self, gene_symbol: str) -> Optional[ChEMBLCompounds]:
        """Return quantitative bioactivity data for compounds active on gene_symbol's target."""
        target_id = await self._resolve_target(gene_symbol)
        if not target_id:
            logger.debug("ChEMBL: no human single-protein target found for %s", gene_symbol)
            return None

        activities = await self._fetch_activities(target_id)
        if not activities:
            return ChEMBLCompounds(
                gene_symbol=gene_symbol,
                target_chembl_id=target_id,
                total_active_compounds=0,
                best_pchembl=None,
                compounds=[],
            )

        best = max(a.pchembl_value for a in activities)
        return ChEMBLCompounds(
            gene_symbol=gene_symbol,
            target_chembl_id=target_id,
            total_active_compounds=len(activities),
            best_pchembl=round(best, 2),
            compounds=activities[:20],
        )

    async def _resolve_target(self, gene_symbol: str) -> Optional[str]:
        """Find the ChEMBL target ID for a human gene symbol."""
        try:
            async with _SEMAPHORE:
                resp = await self._client.get(
                    _TARGET_SEARCH_URL,
                    params={"q": gene_symbol, "format": "json", "limit": 10},
                    timeout=20.0,
                )
                resp.raise_for_status()

            data = resp.json()
            for target in (data.get("targets") or []):
                # Filter: single protein, Homo sapiens
                if (
                    target.get("target_type") == "SINGLE PROTEIN"
                    and target.get("organism") == "Homo sapiens"
                ):
                    return target.get("target_chembl_id")
            return None

        except Exception as exc:
            logger.warning("ChEMBL target search failed for %s: %s", gene_symbol, exc)
            return None

    async def _fetch_activities(self, target_id: str) -> list[ChEMBLActivity]:
        """Fetch IC50/Ki/Kd activities with a pChEMBL value, sorted by potency."""
        try:
            async with _SEMAPHORE:
                resp = await self._client.get(
                    _ACTIVITY_URL,
                    params={
                        "target_chembl_id": target_id,
                        "standard_type__in": "IC50,Ki,Kd,EC50",
                        "pchembl_value__isnull": "false",
                        "limit": 100,
                        "format": "json",
                    },
                    timeout=30.0,
                )
                resp.raise_for_status()

            data = resp.json()
            raw = data.get("activities") or []

            activities: list[ChEMBLActivity] = []
            seen: set[str] = set()

            for a in raw:
                mol_id = a.get("molecule_chembl_id", "")
                pchembl = a.get("pchembl_value")
                if not mol_id or pchembl is None:
                    continue
                try:
                    pchembl_f = float(pchembl)
                except (ValueError, TypeError):
                    continue

                # Deduplicate by molecule: keep best pChEMBL value per compound
                if mol_id in seen:
                    continue
                seen.add(mol_id)

                activities.append(
                    ChEMBLActivity(
                        molecule_chembl_id=mol_id,
                        molecule_name=a.get("molecule_pref_name"),
                        standard_type=a.get("standard_type", ""),
                        pchembl_value=round(pchembl_f, 2),
                        assay_description=a.get("assay_description", "")[:120] if a.get("assay_description") else None,
                    )
                )

            # Sort by potency descending
            activities.sort(key=lambda x: x.pchembl_value, reverse=True)
            return activities

        except Exception as exc:
            logger.warning("ChEMBL activity fetch failed for %s: %s", target_id, exc)
            return []
