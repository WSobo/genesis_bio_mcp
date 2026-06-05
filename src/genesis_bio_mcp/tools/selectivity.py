"""Orchestrate a compound's selectivity / off-target profile.

Combines two honestly-separated layers:
- **Measured** off-targets from ChEMBL (keyless) — what the compound *demonstrably* binds.
- **Predicted** off-targets from the kinome corpus (optional) — kinases hit by chemically-similar
  known corpus binders (SEA-style analogy), which also covers novel structures ChEMBL hasn't seen.

Returns ``None`` only when the SMILES is unparseable; a valid structure with no measured data still
yields a profile (with an explanatory note), so the tool reports "no data" honestly rather than erroring.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from genesis_bio_mcp.corpus import predict_off_targets_by_similarity
from genesis_bio_mcp.models import OffTarget, PredictedOffTarget, SelectivityProfile
from genesis_bio_mcp.tools.cheminformatics import morgan_fp_bits, standardize_structure

if TYPE_CHECKING:
    import asyncpg

    from genesis_bio_mcp.clients.chembl import ChEMBLClient

logger = logging.getLogger(__name__)

# A non-primary target engaged at pChEMBL ≥ this (sub-100 nM) counts as a *potent* off-target.
# An absolute threshold (not a window around the primary) so an exceptionally potent primary —
# e.g. dasatinib's ABL1 at pChEMBL 10 — doesn't mask a wide off-target footprint.
_POTENT_PCHEMBL = 7.0


def _selectivity_label(measured: list[OffTarget]) -> tuple[str, str | None]:
    """Qualitative selectivity from the measured set + a primary-target string."""
    if not measured:
        return "Unknown (no measured data)", None
    primary = measured[0]  # measured is sorted by best pChEMBL desc
    primary_str = (
        f"{primary.gene_symbol or primary.target_pref_name} (pChEMBL {primary.best_pchembl:.1f})"
    )
    n_potent_off = sum(1 for o in measured[1:] if o.best_pchembl >= _POTENT_PCHEMBL)
    label = (
        "Selective"
        if n_potent_off <= 1
        else "Moderately selective"
        if n_potent_off <= 4
        else "Promiscuous"
    )
    return label, primary_str


async def assess_selectivity(
    smiles: str,
    *,
    chembl: ChEMBLClient,
    corpus_pool: asyncpg.Pool | None = None,
) -> SelectivityProfile | None:
    """Build a compound's off-target / selectivity profile. ``None`` iff the SMILES is unparseable."""
    std = standardize_structure(smiles)
    if std is None or not std.inchikey:
        return None

    mol_id, measured = await chembl.get_off_targets(std.inchikey)
    label, primary = _selectivity_label(measured)

    predicted: list[PredictedOffTarget] = []
    corpus_used = corpus_pool is not None
    if corpus_pool is not None:
        fp = morgan_fp_bits(smiles)
        if fp:
            try:
                rows = await predict_off_targets_by_similarity(corpus_pool, fp)
            except Exception as exc:  # never let a corpus issue break the tool
                logger.warning("Corpus off-target prediction failed: %s", exc)
                rows = []
            measured_genes = {o.gene_symbol for o in measured if o.gene_symbol}
            for r in rows:
                gene = r.get("gene_symbol")
                if not gene or gene in measured_genes:
                    continue  # surface only *new* predicted off-targets, not measured ones
                ap = r.get("analog_best_pchembl")
                predicted.append(
                    PredictedOffTarget(
                        gene_symbol=gene,
                        uniprot_accession=r.get("uniprot_accession"),
                        max_tanimoto=round(float(r["max_tanimoto"]), 3),
                        analog_best_pchembl=round(float(ap), 2) if ap is not None else None,
                        n_analogs=int(r["n_analogs"]),
                    )
                )

    return SelectivityProfile(
        query_smiles=std.input_smiles,
        molecule_chembl_id=mol_id,
        primary_target=primary,
        selectivity_label=label,
        measured_off_targets=measured,
        predicted_off_targets=predicted,
        corpus_used=corpus_used,
    )
