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

# ---------------------------------------------------------------------------
# Scoring constants
#
# These are the two tunable parameters most likely to need adjustment as
# ground-truth target assessments accumulate. Both are intentionally modest
# to avoid over-fitting to any single target class.
#
# LINEAGE_MATCH_FACTOR (1.2):
#   Multiplicative bonus applied to the DepMap dependency contribution when
#   the query indication matches a top dependent lineage. Rescues lineage-
#   restricted oncogenes (e.g. BRAF in melanoma) from being penalized for
#   low pan-cancer dependency. Capped so the DepMap axis never exceeds its
#   2.0 max. Range to consider: 1.1–1.5. Above ~1.5 starts to over-reward
#   targets with any incidental lineage overlap.
#
# PUBCHEM_MIN_COMPOUNDS (5):
#   Minimum PubChem active compound count required to contribute a chemical
#   matter score. Below this threshold (and with no ChEMBL potency data),
#   the target is reported as "no usable chemical matter." The threshold
#   filters out single-assay hits and data entry noise. Range to consider:
#   3–10. Below 3 risks counting artifacts; above 10 may be too conservative
#   for early-stage targets.
# ---------------------------------------------------------------------------

LINEAGE_MATCH_FACTOR = 1.2
PUBCHEM_MIN_COMPOUNDS = 5

# OT_CLINICALLY_VALIDATED_FLOOR (3.25):
#   Minimum OT contribution for targets with approved drugs (known_drug_score > 0.9).
#   OT's overall_score is a composite across 6-8 evidence datatypes (genetic, somatic,
#   known_drug, literature, expression, pathways, animal models, etc.). For biologics
#   targets, expression/pathway/animal model datatypes score moderately (~0.4-0.7),
#   pulling the composite below 0.7 even when known_drug and literature are ~0.99.
#   Approved drugs are the gold standard of target-disease validation; a composite
#   score of 0.63-0.64 for TNF/RA and HER2/breast cancer understates that certainty.
#   The floor ensures pharmacological validation is never washed out by weak-but-
#   populated orthogonal evidence types.
#
#   Value exceeds the OT axis natural max (3.0) by design: the excess 0.25 represents
#   certainty from approved drug evidence that OT's multi-datatype average cannot fully
#   express. Calibrated so HER2/breast cancer (OT base 1.902 → 7.65 High) and TNF/RA
#   (OT base 1.926 → 7.10 High) reach the correct tier. Range to consider: 3.0–3.5.
#   Below 3.15 leaves TNF at Medium; above 3.5 risks over-scoring targets with one
#   approved drug and a weak OT overall.
OT_CLINICALLY_VALIDATED_FLOOR = 3.25


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
        (protein, protein_err, t_uniprot),
        (disease_assoc, da_err, t_ot),
        (cancer_dep, cd_err, t_depmap),
        (gwas_ev, gw_err, t_gwas),
        (compounds, co_err, t_pubchem),
        (chembl_compounds, chembl_err, t_chembl),
    ) = await asyncio.gather(
        _safe_timed("uniprot", uniprot.get_protein(symbol)),
        _safe_timed("open_targets", open_targets.get_association(symbol, indication)),
        _safe_timed("depmap", depmap.get_essentiality(symbol)),
        _safe_timed("gwas", gwas.get_evidence(symbol, indication, ncbi_gene_id=ncbi_id)),
        _safe_timed("pubchem", pubchem.get_compounds(symbol)),
        _safe_timed("chembl", chembl.get_compounds(symbol)),
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
        disease_assoc,
        cancer_dep,
        gwas_ev,
        compounds,
        protein,
        chembl_compounds,
        indication=indication,
    )
    priority_tier = _tier(priority_score)

    # Confidence scoring: quantify data completeness and flag proxy sources
    _CORE_SOURCES = ["uniprot", "open_targets", "depmap", "gwas", "pubchem", "chembl"]
    filled = sum(1 for s in _CORE_SOURCES if s not in data_gaps)
    data_coverage_pct = round(filled / len(_CORE_SOURCES) * 100, 1)

    proxy_data_flags: dict[str, bool] = {}
    if cancer_dep is not None and "DepMap Chronos" not in cancer_dep.data_source:
        proxy_data_flags["depmap"] = True
    if chembl_compounds is None and compounds is not None:
        proxy_data_flags["compounds"] = True

    missing_fraction = 1.0 - data_coverage_pct / 100.0
    score_bound = round(priority_score * missing_fraction * 0.5, 2)
    score_confidence_interval: tuple[float, float] | None = None
    if score_bound > 0:
        score_confidence_interval = (
            round(max(0.0, priority_score - score_bound), 2),
            round(min(10.0, priority_score + score_bound), 2),
        )

    evidence_summary = _build_summary(
        symbol,
        indication,
        disease_assoc,
        cancer_dep,
        gwas_ev,
        compounds,
        chembl_compounds,
        ot_error=da_err,
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

    api_latency_s: dict[str, float] = {
        "uniprot": round(t_uniprot, 2),
        "open_targets": round(t_ot, 2),
        "depmap": round(t_depmap, 2),
        "gwas": round(t_gwas, 2),
        "pubchem": round(t_pubchem, 2),
        "chembl": round(t_chembl, 2),
    }
    logger.info(
        "API latencies for %s/%s: %s",
        symbol,
        indication,
        {k: f"{v:.2f}s" for k, v in api_latency_s.items()},
    )

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
        data_coverage_pct=data_coverage_pct,
        proxy_data_flags=proxy_data_flags,
        score_confidence_interval=score_confidence_interval,
        api_latency_s=api_latency_s,
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


async def _safe_timed(name: str, coro, fallback=None) -> tuple[Any, str | None, float]:
    """Like _safe(), but also returns wall-clock seconds for the awaited coroutine."""
    t0 = asyncio.get_running_loop().time()
    result, err = await _safe(coro, fallback)
    return result, err, asyncio.get_running_loop().time() - t0


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
    indication: str = "",
) -> float:
    score = 0.0

    # Open Targets association (max 3.0; floor raised for clinically validated targets)
    # When known_drug_score > 0.9 AND both genetic_association and somatic_mutation
    # subscores are null, floor the OT contribution at OT_CLINICALLY_VALIDATED_FLOOR.
    #
    # The two-gate condition matters: genetic/somatic null means the OT composite is
    # artificially low because those evidence classes are non-applicable for this
    # target mechanism (CNV amplification, cytokine inhibition) — not because the
    # evidence is weak. When somatic or genetic IS populated (e.g. BRAF V600E), the
    # OT composite accurately reflects the full evidence base and should not be floored.
    if disease_assoc:
        ot_contrib = disease_assoc.overall_score * 3.0
        if (
            (disease_assoc.known_drug_score or 0.0) > 0.9
            and disease_assoc.genetic_association_score is None
            and disease_assoc.somatic_mutation_score is None
        ):
            ot_contrib = max(ot_contrib, OT_CLINICALLY_VALIDATED_FLOOR)
        score += ot_contrib

    # DepMap cancer dependency (max 2.0)
    # Apply 0.7x confidence discount when using OT somatic mutation proxy instead of real CRISPR data.
    # Apply 1.2x lineage bonus when the query indication matches a top dependent lineage
    # (e.g. BRAF 9% overall dependency is concentrated in melanoma — penalizing it equally to a
    # pan-cancer 9% target misrepresents the biology for a melanoma query).
    if cancer_dep and not cancer_dep.pan_essential:
        is_real_depmap = "DepMap Chronos" in cancer_dep.data_source
        confidence = 1.0 if is_real_depmap else 0.7
        lineage_match = indication and any(
            indication.lower() in lin.lower() or lin.lower() in indication.lower()
            for lin in (cancer_dep.top_dependent_lineages or [])
        )
        lineage_factor = LINEAGE_MATCH_FACTOR if lineage_match else 1.0
        dep_contribution = min(
            cancer_dep.fraction_dependent_lines * 2.0 * confidence * lineage_factor, 2.0
        )
        score += dep_contribution
    elif cancer_dep and cancer_dep.pan_essential:
        # Pan-essential genes have narrow therapeutic windows — cap contribution
        score += 0.5

    # GWAS evidence (max 2.0)
    # Cap at 3: ≥3 replicated trait hits = full credit. Keeps score stable whether
    # the fetch returned 3 or 30 hits — pagination differences don't affect scoring.
    if gwas_ev:
        score += min(gwas_ev.total_associations, 3) / 3 * 2.0

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
    elif compounds and compounds.total_active_compounds >= PUBCHEM_MIN_COMPOUNDS:
        score += min(compounds.total_active_compounds, 100) / 100 * 1.5
    # PubChem count < PUBCHEM_MIN_COMPOUNDS with no ChEMBL potency data → no usable chemical matter

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
    ot_error: str | None = None,
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
    elif ot_error:
        parts.append(
            f"[OT UNAVAILABLE — score based on ChEMBL/DepMap/GWAS only] "
            f"Open Targets API error for {symbol}/{indication} ({ot_error[:80]}). "
            f"Score may understate target-disease association strength."
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
            lineage_match = indication and any(
                indication.lower() in lin.lower() or lin.lower() in indication.lower()
                for lin in (cancer_dep.top_dependent_lineages or [])
            )
            bonus_note = (
                f" [DepMap score boosted {LINEAGE_MATCH_FACTOR}× — indication matches top lineage]"
                if lineage_match
                else ""
            )
            parts.append(
                f"DepMap CRISPR data show dependency in {pct}% of cancer lines"
                + (f", highest in {top}" if top else "")
                + bonus_note
                + "."
            )

    if gwas_ev and gwas_ev.total_associations > 0:
        p_str = f"{gwas_ev.strongest_p_value:.2e}" if gwas_ev.strongest_p_value else "N/A"
        parts.append(
            f"GWAS Catalog links {gwas_ev.total_associations} variants near {symbol} "
            f"to '{indication}'-related traits (strongest p={p_str})."
        )
        # Causal caveat: high GWAS signal without known-drug evidence AND sparse literature
        # suggests LD noise, comorbidity confounding, or indirect drug effects (e.g. cancer drugs
        # affecting metabolic markers in T2D patient registries).
        # Three-gate condition uses literature_mining_score as the causal discriminator:
        # - FTO/obesity: literature ~0.98 (massive validated biology) → gate suppresses caveat ✓
        # - BRAF/T2D: literature ~0.07 (no real BRAF-T2D biology) → gate fires caveat ✓
        # genetic_association_score was NOT used: OT accumulates GWAS Catalog counts and both
        # targets score ~0.7–0.8, measuring signal volume rather than causal validation.
        # Threshold 0.15: separates replicated functional loci (>0.5) from GWAS noise (<0.1).
        # Tuning range: 0.1–0.25.
        known_drug_score = (disease_assoc.known_drug_score or 0.0) if disease_assoc else 0.0
        literature_score = (disease_assoc.literature_mining_score or 0.0) if disease_assoc else 0.0
        if gwas_ev.total_associations >= 5 and known_drug_score < 0.1 and literature_score < 0.15:
            ot_score = disease_assoc.overall_score if disease_assoc else 0.0
            parts.append(
                f"Caution: {gwas_ev.total_associations} GWAS hits near {symbol} for {indication} "
                f"with no known-drug evidence and sparse literature support "
                f"(OT drug score={known_drug_score:.2f}, literature={literature_score:.2f}, "
                f"overall={ot_score:.2f}). This pattern is consistent with LD with a nearby causal "
                f"variant, comorbidity confounding, or indirect drug effects — not direct target biology. "
                f"Functional validation is recommended before treating GWAS signal as target evidence."
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
    elif compounds and compounds.total_active_compounds >= PUBCHEM_MIN_COMPOUNDS:
        parts.append(
            f"PubChem reports {compounds.total_active_compounds} active compounds against {symbol}, "
            f"indicating {'strong' if compounds.total_active_compounds > 50 else 'emerging'} druggability."
        )
    elif compounds:
        n = compounds.total_active_compounds
        parts.append(
            f"No usable chemical matter found for {symbol} "
            f"(ChEMBL: no potency data; PubChem: {n} hit{'s' if n != 1 else ''} — "
            f"below tractability threshold of {PUBCHEM_MIN_COMPOUNDS} active compounds)."
        )

    return (
        " ".join(parts)
        if parts
        else f"Insufficient data to summarize evidence for {symbol} in {indication}."
    )
