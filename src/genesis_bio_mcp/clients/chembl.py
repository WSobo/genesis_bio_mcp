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

import httpx

from genesis_bio_mcp.config.settings import settings
from genesis_bio_mcp.models import ChEMBLActivity, ChEMBLCompounds

logger = logging.getLogger(__name__)

_BASE = "https://www.ebi.ac.uk/chembl/api/data"
_TARGET_SEARCH_URL = f"{_BASE}/target/search"
_ACTIVITY_URL = f"{_BASE}/activity"
_MECHANISM_URL = f"{_BASE}/mechanism"

_SEMAPHORE = asyncio.Semaphore(settings.chembl_semaphore_limit)


class ChEMBLClient:
    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client

    async def get_compounds(self, gene_symbol: str) -> ChEMBLCompounds | None:
        """Return quantitative bioactivity data for compounds active on gene_symbol's target."""
        target_id = await self._resolve_target(gene_symbol)
        if not target_id:
            logger.debug("ChEMBL: no human single-protein target found for %s", gene_symbol)
            return None

        # Activities (potency) and mechanisms (clinical precedent) are
        # independent endpoints — fetch them concurrently.
        activities, mechanism = await asyncio.gather(
            self._fetch_activities(target_id),
            self._fetch_mechanisms(target_id),
        )
        max_phase, moas, action_types = mechanism
        if not activities:
            return ChEMBLCompounds(
                gene_symbol=gene_symbol,
                target_chembl_id=target_id,
                total_active_compounds=0,
                best_pchembl=None,
                max_clinical_phase=max_phase,
                mechanisms_of_action=moas,
                action_types=action_types,
                compounds=[],
            )

        best = max(a.pchembl_value for a in activities)
        # Split bests by assay type so the scoring axis can prefer functional/
        # cell-based potency over binding-only (which over-credits scaffold
        # binders for undruggable targets like MYC).
        functional = [a.pchembl_value for a in activities if (a.assay_type or "").upper() == "F"]
        binding = [a.pchembl_value for a in activities if (a.assay_type or "").upper() == "B"]
        return ChEMBLCompounds(
            gene_symbol=gene_symbol,
            target_chembl_id=target_id,
            total_active_compounds=len(activities),
            best_pchembl=round(best, 2),
            best_pchembl_functional=round(max(functional), 2) if functional else None,
            best_pchembl_binding=round(max(binding), 2) if binding else None,
            max_clinical_phase=max_phase,
            mechanisms_of_action=moas,
            action_types=action_types,
            compounds=activities[:20],
        )

    async def _resolve_target(self, gene_symbol: str) -> str | None:
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
            for target in data.get("targets") or []:
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

                # Confidence score is a string in some ChEMBL rows and an int
                # in others; coerce defensively, drop on failure.
                conf_raw = a.get("confidence_score")
                try:
                    conf_score = int(conf_raw) if conf_raw is not None else None
                except (ValueError, TypeError):
                    conf_score = None

                activities.append(
                    ChEMBLActivity(
                        molecule_chembl_id=mol_id,
                        molecule_name=a.get("molecule_pref_name"),
                        standard_type=a.get("standard_type", ""),
                        pchembl_value=round(pchembl_f, 2),
                        assay_description=a.get("assay_description", "")[:120]
                        if a.get("assay_description")
                        else None,
                        assay_type=a.get("assay_type") or None,
                        # ChEMBL's /activity rows expose the species via
                        # ``target_organism``; ``assay_organism`` is almost
                        # always null. Fall back to the latter for assays where
                        # the assay system differs from the target species
                        # (e.g. human protein expressed in insect cells).
                        assay_organism=a.get("target_organism") or a.get("assay_organism") or None,
                        assay_cell_type=a.get("assay_cell_type") or None,
                        bao_format=a.get("bao_label") or a.get("bao_format") or None,
                        confidence_score=conf_score,
                    )
                )

            # Sort by potency descending
            activities.sort(key=lambda x: x.pchembl_value, reverse=True)
            return activities

        except Exception as exc:
            logger.warning("ChEMBL activity fetch failed for %s: %s", target_id, exc)
            return []

    async def _fetch_mechanisms(self, target_id: str) -> tuple[float | None, list[str], list[str]]:
        """Fetch mechanism-of-action records for the target.

        Returns ``(max_clinical_phase, mechanisms_of_action, action_types)``.
        The ChEMBL ``/mechanism`` endpoint carries ``max_phase`` inline, so the
        most advanced clinical stage of any drug acting on this target comes
        from a single call (no per-molecule fan-out). Never raises.
        """
        try:
            async with _SEMAPHORE:
                resp = await self._client.get(
                    _MECHANISM_URL,
                    params={
                        "target_chembl_id": target_id,
                        "limit": 100,
                        "format": "json",
                    },
                    timeout=20.0,
                )
                resp.raise_for_status()
            rows = resp.json().get("mechanisms") or []
        except Exception as exc:
            logger.warning("ChEMBL mechanism fetch failed for %s: %s", target_id, exc)
            return None, [], []

        max_phase: float | None = None
        moas: list[str] = []
        action_types: list[str] = []
        for m in rows:
            phase_raw = m.get("max_phase")
            try:
                phase = float(phase_raw) if phase_raw is not None else None
            except (TypeError, ValueError):
                phase = None
            if phase is not None and (max_phase is None or phase > max_phase):
                max_phase = phase
            moa = m.get("mechanism_of_action")
            if isinstance(moa, str) and moa and moa not in moas:
                moas.append(moa)
            at = m.get("action_type")
            if isinstance(at, str) and at and at not in action_types:
                action_types.append(at)
        return max_phase, moas, action_types
