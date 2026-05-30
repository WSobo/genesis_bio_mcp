"""Fetch ChEMBL bioactivities for the corpus's targets (compounds + activities).

Bounded, REST-based ingestion: for each corpus target (UniProt accession) we resolve its
ChEMBL target id, then pull IC50/Ki/Kd/EC50 activities that carry a pChEMBL value. Each
activity row carries the molecule's ``canonical_smiles``, so compounds and activities come
from one query. RDKit (core dep) standardizes each compound + computes its Morgan fingerprint
— no ML model, no ``ingest`` extras needed for this step.

For a *full* kinome build the EBI ChEMBL dump is faster (REST per-target is rate-limited);
this REST path is the accessible default for a partial / demo build (cap with a target limit).
"""

from __future__ import annotations

import logging

import httpx

from genesis_bio_mcp.ingest.load import ActivityRecord, CompoundRecord
from genesis_bio_mcp.tools.cheminformatics import standardized_compound_fields

logger = logging.getLogger(__name__)

_BASE = "https://www.ebi.ac.uk/chembl/api/data"
_TARGET_URL = f"{_BASE}/target"
_ACTIVITY_URL = f"{_BASE}/activity"


def _as_float(value: object) -> float | None:
    try:
        return float(value) if value is not None else None
    except (TypeError, ValueError):
        return None


def _as_int(value: object) -> int | None:
    try:
        return int(value) if value is not None else None
    except (TypeError, ValueError):
        return None


async def resolve_chembl_target(client: httpx.AsyncClient, uniprot_accession: str) -> str | None:
    """Resolve a UniProt accession to a human single-protein ChEMBL target id, or None."""
    try:
        resp = await client.get(
            _TARGET_URL,
            params={
                "target_components__accession": uniprot_accession,
                "format": "json",
                "limit": 20,
            },
        )
        resp.raise_for_status()
    except Exception as exc:
        logger.warning("ChEMBL target lookup failed for %s: %s", uniprot_accession, exc)
        return None
    for target in resp.json().get("targets") or []:
        if (
            target.get("target_type") == "SINGLE PROTEIN"
            and target.get("organism") == "Homo sapiens"
        ):
            return target.get("target_chembl_id")
    return None


def _parse_activity(row: dict, uniprot_accession: str) -> ActivityRecord | None:
    """Map one ChEMBL /activity row to an ActivityRecord, or None if it lacks key fields."""
    act_id = _as_int(row.get("activity_id"))
    mol_id = row.get("molecule_chembl_id")
    tid = row.get("target_chembl_id")
    if act_id is None or not mol_id or not tid:
        return None
    return ActivityRecord(
        activity_id=act_id,
        molecule_chembl_id=mol_id,
        uniprot_accession=uniprot_accession,
        target_chembl_id=tid,
        standard_type=row.get("standard_type"),
        standard_value=_as_float(row.get("standard_value")),
        standard_units=row.get("standard_units"),
        pchembl_value=_as_float(row.get("pchembl_value")),
        assay_chembl_id=row.get("assay_chembl_id"),
        assay_confidence_score=_as_int(row.get("confidence_score")),
        doc_chembl_id=row.get("document_chembl_id"),
    )


async def fetch_activities_for_target(
    client: httpx.AsyncClient,
    target_chembl_id: str,
    uniprot_accession: str,
    *,
    page_size: int = 1000,
    max_pages: int = 5,
) -> tuple[list[ActivityRecord], dict[str, str | None]]:
    """Fetch IC50/Ki/Kd/EC50 activities (with a pChEMBL value) for a ChEMBL target.

    Follows ChEMBL ``page_meta.next`` pagination up to ``max_pages`` (a generous cap for a
    demo build). Returns ``(activities, smiles_by_molecule)`` — the SMILES map is used to
    build the compound rows (rows missing key fields are dropped).
    """
    records: list[ActivityRecord] = []
    smiles_by_molecule: dict[str, str | None] = {}
    params: dict | None = {
        "target_chembl_id": target_chembl_id,
        "standard_type__in": "IC50,Ki,Kd,EC50",
        "pchembl_value__isnull": "false",
        "format": "json",
        "limit": page_size,
    }
    next_path: str | None = None
    for _ in range(max_pages):
        url = f"{_BASE.split('/chembl')[0]}{next_path}" if next_path else _ACTIVITY_URL
        try:
            resp = await client.get(url, params=None if next_path else params)
            resp.raise_for_status()
        except Exception as exc:
            logger.warning("ChEMBL activity fetch failed for %s: %s", target_chembl_id, exc)
            break
        body = resp.json()
        for row in body.get("activities") or []:
            rec = _parse_activity(row, uniprot_accession)
            if rec is not None:
                records.append(rec)
                smiles_by_molecule.setdefault(rec.molecule_chembl_id, row.get("canonical_smiles"))
        next_path = (body.get("page_meta") or {}).get("next")
        if not next_path:
            break
    return records, smiles_by_molecule


def build_compound_record(molecule_chembl_id: str, smiles: str | None) -> CompoundRecord:
    """Build a CompoundRecord from a SMILES (RDKit canonical + InChIKey + MW + Morgan FP)."""
    fields = standardized_compound_fields(smiles) if smiles else None
    if fields is None:
        return CompoundRecord(molecule_chembl_id, smiles, None, None, None, None)
    canonical, inchi, inchikey, mw, fp_bits = fields
    return CompoundRecord(molecule_chembl_id, canonical, inchi, inchikey, mw, fp_bits)
