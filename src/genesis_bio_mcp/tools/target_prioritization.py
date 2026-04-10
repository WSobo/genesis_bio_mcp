"""Orchestration tool: chains all database clients into a target assessment report."""

from __future__ import annotations

import asyncio
import logging
from typing import TYPE_CHECKING, Any

from genesis_bio_mcp.models import (
    CancerDependency,
    ChEMBLCompounds,
    Compounds,
    DrugHistory,
    GeneResolution,
    GwasEvidence,
    PathwayContext,
    ProteinInfo,
    ProteinInteractome,
    ProteinStructure,
    TargetDiseaseAssociation,
    TargetPrioritizationReport,
)
from genesis_bio_mcp.tools.gene_resolver import resolve_gene

if TYPE_CHECKING:
    from genesis_bio_mcp.clients.alphafold import AlphaFoldClient
    from genesis_bio_mcp.clients.chembl import ChEMBLClient
    from genesis_bio_mcp.clients.clinical_trials import ClinicalTrialsClient
    from genesis_bio_mcp.clients.depmap import DepMapClient
    from genesis_bio_mcp.clients.dgidb import DGIdbClient
    from genesis_bio_mcp.clients.gwas import GwasClient
    from genesis_bio_mcp.clients.open_targets import OpenTargetsClient
    from genesis_bio_mcp.clients.pubchem import PubChemClient
    from genesis_bio_mcp.clients.reactome import ReactomeClient
    from genesis_bio_mcp.clients.string_db import StringDbClient
    from genesis_bio_mcp.clients.uniprot import UniProtClient

logger = logging.getLogger(__name__)


async def prioritize_target(
    gene_symbol: str,
    indication: str,
    *,
    uniprot: UniProtClient,
    open_targets: OpenTargetsClient,
    depmap: DepMapClient,
    gwas: GwasClient,
    pubchem: PubChemClient,
    chembl: ChEMBLClient,
    alphafold: AlphaFoldClient | None = None,
    string_db: StringDbClient | None = None,
    dgidb: DGIdbClient | None = None,
    clinical_trials: ClinicalTrialsClient | None = None,
    reactome: ReactomeClient | None = None,
) -> TargetPrioritizationReport:
    """Run all database queries in parallel and synthesize a target prioritization report.

    Never raises — failed sub-queries are captured in data_gaps/errors.
    """
    symbol = gene_symbol.strip().upper()

    # Resolve gene first (needed for NCBI gene ID in GWAS lookup)
    resolution, resolution_err = await _safe(
        resolve_gene(symbol, uniprot_client=uniprot),
        fallback=GeneResolution(hgnc_symbol=symbol, source="input"),
    )

    ncbi_id = resolution.ncbi_gene_id if resolution else None

    # Use canonical symbol for all downstream lookups (handles HER2→ERBB2, p53→TP53, COX2→PTGS2)
    if resolution and resolution.hgnc_symbol and resolution.hgnc_symbol != symbol:
        logger.info(
            "Resolved alias %s → %s; using canonical symbol for all lookups",
            symbol,
            resolution.hgnc_symbol,
        )
        symbol = resolution.hgnc_symbol

    # Run all independent lookups concurrently
    (
        (protein, protein_err),
        (disease_assoc, da_err),
        (cancer_dep, cd_err),
        (gwas_ev, gw_err),
        (compounds, co_err),
        (chembl_compounds, chembl_err),
    ) = await asyncio.gather(
        _safe(uniprot.get_protein(symbol)),
        _safe(open_targets.get_association(symbol, indication)),
        _safe(depmap.get_essentiality(symbol)),
        _safe(gwas.get_evidence(symbol, indication, ncbi_gene_id=ncbi_id)),
        _safe(pubchem.get_compounds(symbol)),
        _safe(chembl.get_compounds(symbol)),
    )

    data_gaps: list[str] = []
    errors: dict[str, str] = {}

    def _check(name: str, value: Any, err: str | None) -> None:
        if err:
            errors[name] = err
            data_gaps.append(name)
        elif value is None:
            data_gaps.append(name)

    if resolution_err:
        errors["gene_resolution"] = resolution_err

    _check("uniprot", protein, protein_err)
    _check("open_targets", disease_assoc, da_err)
    _check("depmap", cancer_dep, cd_err)
    _check("gwas", gwas_ev, gw_err)
    _check("pubchem", compounds, co_err)
    _check("chembl", chembl_compounds, chembl_err)

    priority_score = _compute_score(
        disease_assoc, cancer_dep, gwas_ev, compounds, protein, chembl_compounds
    )
    priority_tier = _tier(priority_score)
    evidence_summary = _build_summary(
        symbol,
        indication,
        disease_assoc,
        cancer_dep,
        gwas_ev,
        compounds,
        chembl_compounds,
    )

    # Extended mode: structure, interactome, drug history, pathway context
    ext_structure: ProteinStructure | None = None
    ext_interactome: ProteinInteractome | None = None
    ext_drug_history: DrugHistory | None = None
    ext_pathway: PathwayContext | None = None

    if any(c is not None for c in (alphafold, string_db, dgidb, clinical_trials, reactome)):
        uniprot_accession = protein.uniprot_accession if protein else None

        ext_results = await asyncio.gather(
            _safe(alphafold.get_structure(symbol, uniprot_accession=uniprot_accession))
            if alphafold
            else _safe_none(),
            _safe(string_db.get_interactome(symbol)) if string_db else _safe_none(),
            _safe(_fetch_drug_history(symbol, dgidb, clinical_trials))
            if (dgidb or clinical_trials)
            else _safe_none(),
            _safe(reactome.get_pathway_context(symbol)) if reactome else _safe_none(),
        )
        (
            (ext_structure, _),
            (ext_interactome, _),
            (ext_drug_history, _),
            (ext_pathway, _),
        ) = ext_results

    return TargetPrioritizationReport(
        gene_symbol=symbol,
        indication=indication,
        resolution=resolution,
        protein_info=protein,
        disease_association=disease_assoc,
        cancer_dependency=cancer_dep,
        gwas_evidence=gwas_ev,
        compounds=compounds,
        chembl_compounds=chembl_compounds,
        priority_score=round(priority_score, 2),
        priority_tier=priority_tier,
        evidence_summary=evidence_summary,
        data_gaps=data_gaps,
        errors=errors,
        protein_structure=ext_structure,
        protein_interactome=ext_interactome,
        drug_history=ext_drug_history,
        pathway_context=ext_pathway,
    )


async def _safe(coro, fallback=None) -> tuple[Any, str | None]:
    """Await a coroutine, returning (result, error_str). Never raises."""
    try:
        result = await coro
        return result, None
    except Exception as exc:
        logger.warning("Tool sub-query failed: %s", exc)
        return fallback, str(exc)


async def _safe_none() -> tuple[None, None]:
    """Placeholder coroutine for skipped extended-mode lookups."""
    return None, None


async def _fetch_drug_history(
    gene_symbol: str,
    dgidb: Any,
    clinical_trials: Any,
) -> DrugHistory | None:
    """Combine DGIdb + ClinicalTrials results into a DrugHistory object."""
    coros = []
    if dgidb is not None:
        coros.append(dgidb.get_drug_interactions(gene_symbol))
    else:
        coros.append(_return_empty_list())
    if clinical_trials is not None:
        coros.append(clinical_trials.get_trials(gene_symbol))
    else:
        coros.append(_return_empty_tuple())

    drugs, ct_result = await asyncio.gather(*coros)
    ct_trials, ct_counts = ct_result if isinstance(ct_result, tuple) else ([], {})
    approved_count = sum(1 for d in drugs if d.approved)
    return DrugHistory(
        gene_symbol=gene_symbol,
        known_drugs=drugs,
        approved_drug_count=approved_count,
        trial_counts_by_phase=ct_counts,
        recent_trials=ct_trials[:10],
    )


async def _return_empty_list() -> list:
    return []


async def _return_empty_tuple() -> tuple:
    return [], {}


def _compute_score(
    disease_assoc: TargetDiseaseAssociation | None,
    cancer_dep: CancerDependency | None,
    gwas_ev: GwasEvidence | None,
    compounds: Compounds | None,
    protein: ProteinInfo | None,
    chembl_compounds: ChEMBLCompounds | None = None,
) -> float:
    score = 0.0

    # Open Targets association (max 3.0)
    if disease_assoc:
        score += disease_assoc.overall_score * 3.0

    # DepMap cancer dependency (max 2.0)
    # Apply 0.7x confidence discount when using OT somatic mutation proxy instead of real CRISPR data
    if cancer_dep and not cancer_dep.pan_essential:
        is_real_depmap = "DepMap Chronos" in cancer_dep.data_source
        confidence = 1.0 if is_real_depmap else 0.7
        score += cancer_dep.fraction_dependent_lines * 2.0 * confidence
    elif cancer_dep and cancer_dep.pan_essential:
        # Pan-essential genes have narrow therapeutic windows — cap contribution
        score += 0.5

    # GWAS evidence (max 2.0)
    if gwas_ev:
        score += min(gwas_ev.total_associations, 10) / 10 * 2.0

    # Clinical / known-drug evidence (max 1.5)
    # Distinguishes targets with approved drugs from literature-only at the same OT overall_score
    if disease_assoc and disease_assoc.known_drug_score:
        score += disease_assoc.known_drug_score * 1.5

    # Chemical matter (max 1.5)
    # ChEMBL potency-based scoring takes precedence over PubChem count-based scoring
    if chembl_compounds and chembl_compounds.best_pchembl is not None:
        bp = chembl_compounds.best_pchembl
        if bp >= 9:
            score += 1.5  # IC50 ≤ 1 nM — clinical-grade potency
        elif bp >= 7:
            score += 1.0  # IC50 ≤ 100 nM — lead quality
        elif bp >= 5:
            score += 0.5  # IC50 ≤ 10 µM — hit quality
        else:
            score += 0.25  # active but weak
    elif compounds:
        score += min(compounds.total_active_compounds, 100) / 100 * 1.5

    # Protein quality (max 1.5)
    if protein:
        if protein.reviewed:
            score += 0.5
        score += min(len(protein.known_variants), 2) / 2 * 1.0

    return min(score, 10.0)


def _tier(score: float) -> str:
    if score >= 7.0:
        return "High"
    if score >= 4.0:
        return "Medium"
    return "Low"


def _build_summary(
    symbol: str,
    indication: str,
    disease_assoc: TargetDiseaseAssociation | None,
    cancer_dep: CancerDependency | None,
    gwas_ev: GwasEvidence | None,
    compounds: Compounds | None,
    chembl_compounds: ChEMBLCompounds | None = None,
) -> str:
    parts: list[str] = []

    if disease_assoc:
        strength = "strong" if disease_assoc.overall_score >= 0.5 else "modest"
        parts.append(
            f"{symbol} shows {strength} Open Targets association with {indication} "
            f"(score: {disease_assoc.overall_score:.2f}, n={disease_assoc.evidence_count} evidence items)."
        )
        if disease_assoc.known_drug_score and disease_assoc.known_drug_score >= 0.5:
            parts.append(
                f"Open Targets reports strong known-drug evidence for {symbol} "
                f"(score: {disease_assoc.known_drug_score:.2f}), suggesting existing approved or "
                f"clinical-stage therapeutics — likely biologics if small-molecule data is sparse."
            )
    else:
        parts.append(f"No Open Targets association data found for {symbol} in {indication}.")

    if cancer_dep:
        pct = int(cancer_dep.fraction_dependent_lines * 100)
        if cancer_dep.pan_essential:
            parts.append(
                f"CRISPR screens show {symbol} is pan-essential ({pct}% of lines dependent), "
                f"suggesting a potentially narrow therapeutic window."
            )
        else:
            top = ", ".join(cancer_dep.top_dependent_lineages[:3])
            parts.append(
                f"DepMap CRISPR data show dependency in {pct}% of cancer lines"
                + (f", highest in {top}" if top else "")
                + "."
            )

    if gwas_ev and gwas_ev.total_associations > 0:
        p_str = f"{gwas_ev.strongest_p_value:.2e}" if gwas_ev.strongest_p_value else "N/A"
        parts.append(
            f"GWAS Catalog links {gwas_ev.total_associations} variants near {symbol} "
            f"to '{indication}'-related traits (strongest p={p_str})."
        )

    if chembl_compounds and chembl_compounds.best_pchembl is not None:
        bp = chembl_compounds.best_pchembl
        ic50_nm = 10 ** (9 - bp)
        ic50_str = f"{ic50_nm:.1f} nM" if ic50_nm < 1000 else f"{ic50_nm / 1000:.1f} µM"
        potency_label = (
            "clinical-grade" if bp >= 9 else "lead-quality" if bp >= 7 else "hit-quality"
        )
        parts.append(
            f"ChEMBL reports {chembl_compounds.total_active_compounds} compounds with potency data "
            f"against {symbol}; best IC50 ≈ {ic50_str} ({potency_label}, pChEMBL={bp:.1f})."
        )
    elif compounds:
        parts.append(
            f"PubChem reports {compounds.total_active_compounds} active compounds against {symbol}, "
            f"indicating {'strong' if compounds.total_active_compounds > 50 else 'emerging'} druggability."
        )

    return (
        " ".join(parts)
        if parts
        else f"Insufficient data to summarize evidence for {symbol} in {indication}."
    )
