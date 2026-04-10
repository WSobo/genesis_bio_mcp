"""ClinicalTrials.gov API v2 client."""

from __future__ import annotations

import logging

import httpx

from genesis_bio_mcp.models import ClinicalTrial

logger = logging.getLogger(__name__)

_CT_URL = "https://clinicaltrials.gov/api/v2/studies"

_PHASE_NORMALIZE = {
    "PHASE1": "Phase 1",
    "PHASE2": "Phase 2",
    "PHASE3": "Phase 3",
    "PHASE4": "Phase 4",
    "EARLY_PHASE1": "Phase 1 (Early)",
    "NA": "N/A",
}


class ClinicalTrialsClient:
    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client

    async def get_trials(
        self, gene_symbol: str, max_results: int = 50
    ) -> tuple[list[ClinicalTrial], dict[str, int]]:
        """Return trials mentioning a gene target and counts by phase."""
        try:
            resp = await self._client.get(
                _CT_URL,
                params={
                    "query.term": gene_symbol,
                    "query.intr": gene_symbol,
                    "pageSize": max_results,
                    "format": "json",
                    "fields": "NCTId,BriefTitle,Phase,OverallStatus,Condition",
                },
                timeout=20.0,
            )
            resp.raise_for_status()
            data = resp.json()
            return _parse_trials(data)
        except Exception as exc:
            logger.warning("ClinicalTrials.gov request failed for %s: %s", gene_symbol, exc)
            return [], {}


def _parse_trials(data: dict) -> tuple[list[ClinicalTrial], dict[str, int]]:
    studies = data.get("studies", [])
    trials: list[ClinicalTrial] = []
    phase_counts: dict[str, int] = {}

    for study in studies:
        proto = study.get("protocolSection", {})
        id_mod = proto.get("identificationModule", {})
        status_mod = proto.get("statusModule", {})
        design_mod = proto.get("designModule", {})
        conditions_mod = proto.get("conditionsModule", {})

        nct_id = id_mod.get("nctId", "")
        title = id_mod.get("briefTitle", "")
        status = status_mod.get("overallStatus", "UNKNOWN")

        # Phase: design module has phases list
        raw_phases = design_mod.get("phases", []) or []
        if raw_phases:
            raw_phase = raw_phases[0]
            phase = _PHASE_NORMALIZE.get(raw_phase, raw_phase.replace("_", " ").title())
        else:
            phase = "N/A"

        # Primary condition
        conditions = conditions_mod.get("conditions", [])
        indication = conditions[0] if conditions else None

        if nct_id:
            trials.append(
                ClinicalTrial(
                    nct_id=nct_id,
                    title=title[:120],
                    phase=phase,
                    status=status,
                    indication=indication,
                )
            )
            phase_counts[phase] = phase_counts.get(phase, 0) + 1

    return trials, phase_counts
