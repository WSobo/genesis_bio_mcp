"""Orchestration tool: chains all database clients into a target assessment report."""

from __future__ import annotations

import asyncio
import logging
from typing import TYPE_CHECKING, Any

from genesis_bio_mcp.clients.dgidb import _DIRECT_TYPES
from genesis_bio_mcp.config.settings import settings
from genesis_bio_mcp.models import (
    _EXPRESSION_RATING,
    CancerDependency,
    ChEMBLCompounds,
    Compounds,
    DrugHistory,
    DrugInteraction,
    EvidenceAxis,
    GeneResolution,
    GwasEvidence,
    PathwayContext,
    ProteinAtlasReport,
    ProteinInfo,
    ProteinInteractome,
    ProteinStructure,
    TargetComparisonRow,
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
    from genesis_bio_mcp.clients.hpa import HPAClient
    from genesis_bio_mcp.clients.open_targets import OpenTargetsClient
    from genesis_bio_mcp.clients.openfda import OpenFDAClient
    from genesis_bio_mcp.clients.pubchem import PubChemClient
    from genesis_bio_mcp.clients.reactome import ReactomeClient
    from genesis_bio_mcp.clients.string_db import StringDbClient
    from genesis_bio_mcp.clients.uniprot import UniProtClient

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Classification thresholds for the evidence profile
#
# These are *classification boundaries* used to label one axis at a time
# (strong / moderate / weak / none / n-a) — NOT weights, multipliers, or a
# composite score. None of them is tuned to land a particular named target in
# a particular tier; each is a conventional cutoff for its single metric.
# ---------------------------------------------------------------------------

# PUBCHEM_MIN_COMPOUNDS (5):
#   Minimum PubChem active-compound count for the chemical-tractability axis to
#   register any signal when no ChEMBL potency data exists. Below this (and with
#   no ChEMBL data) the axis is "none" — the count filters single-assay hits and
#   data-entry noise from a real chemical-matter signal.
PUBCHEM_MIN_COMPOUNDS = 5

# OT_GENETIC_MONOGENIC_THRESHOLD (0.7):
#   Open Targets ``genetic_association_score`` above which the genetic-evidence
#   axis is rated "strong" even when GWAS Catalog returns no hits. A score >0.7
#   is the signature of monogenic / Mendelian disease (CFTR/CF≈0.99,
#   HBB/sickle-cell≈0.96, HTT/Huntington's≈0.94) — studied via Mendelian linkage
#   and ClinVar, not population GWAS, so an empty GWAS Catalog is expected rather
#   than weak evidence. Polygenic targets (FTO/obesity≈0.55, PCSK9/CHD≈0.65) sit
#   below it and are rated from their actual GWAS hits.
OT_GENETIC_MONOGENIC_THRESHOLD = 0.7

# _GWAS_TRAIT_RELEVANCE_MIN (0.15):
#   Minimum token-set Jaccard between the indication and a GWAS hit's trait label
#   for that hit to count as *trait-relevant*. This is a relevance/correctness
#   filter, not a score weight: it stops a hit that merely substring-matched
#   something incidental in the Catalog (e.g. a TP53-locus sex-hormone-binding-
#   globulin study for a Li-Fraumeni query) from being counted as genetic
#   evidence for the indication. It is a coarse string-overlap heuristic — it
#   classifies hits as on/off-trait, nothing more.
_GWAS_TRAIT_RELEVANCE_MIN = 0.15


def _max_trait_relevance(associations: list, indication: str) -> float:
    """Highest token-set Jaccard between the indication and any hit's trait label."""
    if not indication or not associations:
        return 0.0
    ind_tokens = _tokenize_for_relevance(indication)
    if not ind_tokens:
        return 0.0
    best = 0.0
    for hit in associations:
        trait = getattr(hit, "trait", None) or ""
        h_tokens = _tokenize_for_relevance(trait)
        if not h_tokens:
            continue
        intersection = ind_tokens & h_tokens
        union = ind_tokens | h_tokens
        if union:
            best = max(best, len(intersection) / len(union))
    return best


def _trait_relevant_gwas_count(gwas_ev: GwasEvidence | None, indication: str) -> int:
    """GWAS hit count that actually matches the indication.

    Single source of truth shared by the per-axis genetic rating
    (``_build_axes``) and the prose summary (``_build_summary``) so the two can
    never disagree. Returns 0 when the GWAS lookup fell back to a non-exact
    trait query, or when the returned hits' trait labels are off-indication
    (token-set Jaccard below ``_GWAS_TRAIT_RELEVANCE_MIN``).
    """
    if gwas_ev is None or gwas_ev.total_associations <= 0:
        return 0
    if "no exact-trait match" in (gwas_ev.trait_query or ""):
        return 0
    if (
        bool(gwas_ev.associations)
        and _max_trait_relevance(gwas_ev.associations, indication) < _GWAS_TRAIT_RELEVANCE_MIN
    ):
        return 0
    return gwas_ev.total_associations


# Stop-word set for trait-relevance tokenization. These tokens appear in too
# many disease/trait labels to discriminate (e.g. "disease", "syndrome") and
# would inflate Jaccard scores spuriously if kept.
_TRAIT_STOPWORDS: frozenset[str] = frozenset(
    {"disease", "syndrome", "disorder", "the", "of", "and", "a", "an", "in", "for"}
)


def _tokenize_for_relevance(text: str) -> set[str]:
    """Lowercase, split on non-alphanumerics, drop stopwords + length-1 tokens."""
    import re

    raw = re.split(r"[^a-z0-9]+", text.lower())
    return {t for t in raw if len(t) > 1 and t not in _TRAIT_STOPWORDS}


def build_comparison_row(
    gene_symbol: str,
    report: TargetPrioritizationReport | BaseException,
) -> TargetComparisonRow:
    """Build one comparison-table row from a prioritization report.

    Shared by ``compare_targets`` (server) and the workflow agent so the two
    surfaces stay in lockstep. Accepts either a finished report or the
    ``Exception`` returned by ``asyncio.gather(..., return_exceptions=True)``.

    PubChem (chemical-space breadth) and ChEMBL (peak measured potency) are kept
    as distinct fields rather than collapsed into a single "compounds" number.
    """
    if isinstance(report, BaseException):
        return TargetComparisonRow(
            gene_symbol=gene_symbol.upper(),
            evidence_profile=[],
            data_gaps=["all"],
            evidence_summary=f"Query failed: {report}",
        )

    cd = report.cancer_dependency
    depmap_pct = int(cd.fraction_dependent_lines * 100) if cd else None
    depmap_real = cd is not None and "DepMap Chronos" in cd.data_source
    chembl = report.chembl_compounds
    return TargetComparisonRow(
        gene_symbol=report.gene_symbol,
        evidence_profile=report.evidence_profile,
        ot_score=report.disease_association.overall_score if report.disease_association else None,
        depmap_pct=depmap_pct,
        depmap_real_data=depmap_real,
        compound_count=report.compounds.total_active_compounds if report.compounds else None,
        chembl_count=chembl.total_active_compounds if chembl else None,
        best_pchembl=chembl.best_pchembl if chembl else None,
        gwas_count=report.gwas_evidence.total_associations if report.gwas_evidence else None,
        data_gaps=report.data_gaps,
        evidence_summary=report.evidence_summary,
    )


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
    openfda: OpenFDAClient | None = None,
    reactome: ReactomeClient | None = None,
    hpa: HPAClient | None = None,
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

    # Extended mode: structure, interactome, drug history, pathway context,
    # tissue expression (HPA). All are independent so they fan out in parallel.
    # Must happen before scoring so the HPA report can feed the expression axis.
    ext_structure: ProteinStructure | None = None
    ext_interactome: ProteinInteractome | None = None
    ext_drug_history: DrugHistory | None = None
    ext_pathway: PathwayContext | None = None
    ext_protein_atlas: ProteinAtlasReport | None = None

    if any(c is not None for c in (alphafold, string_db, dgidb, clinical_trials, reactome, hpa)):
        uniprot_accession = protein.uniprot_accession if protein else None

        ext_results = await asyncio.gather(
            _safe(alphafold.get_structure(symbol, uniprot_accession=uniprot_accession))
            if alphafold
            else _safe_none(),
            _safe(string_db.get_interactome(symbol)) if string_db else _safe_none(),
            _safe(_fetch_drug_history(symbol, dgidb, clinical_trials, openfda=openfda))
            if (dgidb or clinical_trials)
            else _safe_none(),
            _safe(reactome.get_pathway_context(symbol)) if reactome else _safe_none(),
            _safe(hpa.get_report(symbol)) if hpa else _safe_none(),
        )
        (
            (ext_structure, _),
            (ext_interactome, _),
            (ext_drug_history, _),
            (ext_pathway, _),
            (ext_protein_atlas, _),
        ) = ext_results

    evidence_profile = _build_evidence_profile(
        disease_assoc,
        cancer_dep,
        gwas_ev,
        compounds,
        protein,
        chembl_compounds,
        indication=indication,
        protein_atlas=ext_protein_atlas,
    )

    # Data completeness + proxy flags (kept — real provenance, no composite score)
    _CORE_SOURCES = ["uniprot", "open_targets", "depmap", "gwas", "pubchem", "chembl"]
    filled = sum(1 for s in _CORE_SOURCES if s not in data_gaps)
    data_coverage_pct = round(filled / len(_CORE_SOURCES) * 100, 1)

    proxy_data_flags: dict[str, bool] = {}
    if cancer_dep is not None and "DepMap Chronos" not in cancer_dep.data_source:
        proxy_data_flags["depmap"] = True
    if chembl_compounds is None and compounds is not None:
        proxy_data_flags["compounds"] = True

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
        evidence_profile=evidence_profile,
        evidence_summary=evidence_summary,
        data_gaps=data_gaps,
        errors=errors,
        data_coverage_pct=data_coverage_pct,
        proxy_data_flags=proxy_data_flags,
        api_latency_s=api_latency_s,
        protein_structure=ext_structure,
        protein_interactome=ext_interactome,
        drug_history=ext_drug_history,
        pathway_context=ext_pathway,
        protein_atlas=ext_protein_atlas,
    )


async def _safe(coro, fallback=None, *, timeout: float | None = None) -> tuple[Any, str | None]:
    """Await a coroutine, returning (result, error_str). Never raises.

    ``timeout`` (seconds, default ``settings.prioritize_source_budget_secs``) bounds
    the await via ``asyncio.wait_for`` so a single slow/hung sub-query can't stall the
    whole assessment. On expiry the coroutine is cancelled — releasing any client
    semaphore and HTTP connection it held — and we degrade to ``fallback`` with a
    "timed out" marker that surfaces as a data gap.
    """
    if timeout is None:
        timeout = settings.prioritize_source_budget_secs
    try:
        result = await asyncio.wait_for(coro, timeout)
        return result, None
    except TimeoutError:
        logger.warning("Tool sub-query timed out after %.0fs", timeout)
        return fallback, f"timed out after {timeout:.0f}s"
    except Exception as exc:
        logger.warning("Tool sub-query failed: %s", exc)
        return fallback, str(exc)


async def _safe_timed(
    name: str, coro, fallback=None, *, timeout: float | None = None
) -> tuple[Any, str | None, float]:
    """Like _safe(), but also returns wall-clock seconds for the awaited coroutine."""
    t0 = asyncio.get_running_loop().time()
    result, err = await _safe(coro, fallback, timeout=timeout)
    return result, err, asyncio.get_running_loop().time() - t0


async def _safe_none() -> tuple[None, None]:
    """Placeholder coroutine for skipped extended-mode lookups."""
    return None, None


async def _fetch_drug_history(
    gene_symbol: str,
    dgidb: Any,
    clinical_trials: Any,
    openfda: Any = None,
) -> DrugHistory | None:
    """Combine DGIdb + ClinicalTrials + OpenFDA results into a DrugHistory object."""
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
    drugs = await attach_safety_signals(drugs, openfda=openfda)
    approved_count = sum(1 for d in drugs if d.approved)
    return DrugHistory(
        gene_symbol=gene_symbol,
        known_drugs=drugs,
        approved_drug_count=approved_count,
        trial_counts_by_phase=ct_counts,
        recent_trials=ct_trials[:10],
    )


# Cap safety lookups at the top-N approved drugs so adding OpenFDA to
# get_drug_history doesn't blow latency budget on targets with 30+ approved
# drugs (e.g. TNF, EGFR). Five is the display cap; we query a wider pool
# so we can re-rank by FAERS report volume after lookup.
_SAFETY_LOOKUP_CAP = 5
_SAFETY_LOOKUP_POOL = 10


async def attach_safety_signals(
    drugs: list[DrugInteraction],
    *,
    openfda: OpenFDAClient | None,
) -> list[DrugInteraction]:
    """Populate ``.safety`` on the top-N approved direct-engagement drugs.

    Non-approved drugs are left untouched — FAERS data for investigational
    compounds is sparse and often misleading. Failures are swallowed; a drug
    whose OpenFDA lookup errors simply carries ``safety=None``.

    Selection (Bugs F.1 + F.2 — v0.3.3):
    - Strict direct-engagement filter: only drugs whose DGIdb interaction_type
      is in ``_DIRECT_TYPES`` are eligible. The v0.3.2 fallback for untyped
      drugs leaked non-target co-mentions (atezolizumab on EGFR, ACYCLOVIR
      on ABL1) into the panel. DGIdb leaves ``interaction_type=None`` for
      exactly those spurious associations, so allowing untyped drugs in
      defeats the filter.
    - Pre-rank candidates by (phase desc, source-count desc, name asc).
      Source count proxies for "real curated interactor" — more DGIdb
      databases attesting to the interaction = higher confidence.
    - Query a wider pool (``_SAFETY_LOOKUP_POOL``=10) than we display, then
      re-rank by ``total_reports`` so well-attested drugs (decades of FAERS
      reports — IMATINIB) surface ahead of recent approvals with sparse
      FAERS history (ASCIMINIB) that previously won on alphabetical sort.
    """
    if openfda is None or not drugs:
        return drugs

    candidates = [
        d
        for d in drugs
        if d.approved and d.interaction_type and d.interaction_type.lower() in _DIRECT_TYPES
    ]
    if not candidates:
        return drugs

    candidates.sort(
        key=lambda d: (-(d.phase or 4), -len(d.sources), d.drug_name.lower()),
    )
    pool = candidates[:_SAFETY_LOOKUP_POOL]

    signals = await asyncio.gather(
        *(openfda.get_safety_signals(d.drug_name) for d in pool),
        return_exceptions=True,
    )

    drug_signal_pairs: list[tuple[DrugInteraction, Any]] = []
    for drug, signal in zip(pool, signals, strict=True):
        if isinstance(signal, Exception):
            logger.warning("OpenFDA lookup raised for %s: %s", drug.drug_name, signal)
            continue
        if signal is None:
            continue
        drug_signal_pairs.append((drug, signal))

    # Re-rank by FAERS report volume — established drugs with thousands of
    # reports rank above newcomers with double-digit reports.
    drug_signal_pairs.sort(key=lambda ds: ds[1].total_reports or 0, reverse=True)
    selected = drug_signal_pairs[:_SAFETY_LOOKUP_CAP]

    if not selected:
        return drugs

    safety_by_name: dict[str, Any] = {drug.drug_name.lower(): sig for drug, sig in selected}
    out: list[DrugInteraction] = []
    for d in drugs:
        sig = safety_by_name.get(d.drug_name.lower())
        if sig is None:
            out.append(d)
        else:
            out.append(d.model_copy(update={"safety": sig}))
    return out


async def _return_empty_list() -> list:
    return []


async def _return_empty_tuple() -> tuple:
    return [], {}


def _build_evidence_profile(
    disease_assoc: TargetDiseaseAssociation | None,
    cancer_dep: CancerDependency | None,
    gwas_ev: GwasEvidence | None,
    compounds: Compounds | None,
    protein: ProteinInfo | None,
    chembl_compounds: ChEMBLCompounds | None = None,
    indication: str = "",
    protein_atlas: ProteinAtlasReport | None = None,
) -> list[EvidenceAxis]:
    """Classify each evidence axis independently into strong / moderate / weak / none / n-a.

    There is deliberately **no** composite score: each axis is rated only from its
    own source metric, using documented conventional boundaries — never weighted,
    multiplied, or summed. ``none`` means data was present but showed no signal;
    ``n/a`` means the source returned nothing (a gap). Keeping the axes separate
    preserves the cross-axis trade-off (genetics vs tractability vs chemistry are
    distinct go/no-go gates) that a single number would erase.
    """
    axes: list[EvidenceAxis] = []

    # --- Disease association (Open Targets overall_score) --------------------
    if disease_assoc is not None:
        s = disease_assoc.overall_score
        rating = "strong" if s >= 0.5 else "moderate" if s >= 0.2 else "weak" if s > 0 else "none"
        axes.append(
            EvidenceAxis(
                name="Disease association",
                rating=rating,
                detail=f"OT overall {s:.2f}, n={disease_assoc.evidence_count} evidence items",
                value=s,
            )
        )
    else:
        axes.append(
            EvidenceAxis(name="Disease association", rating="n/a", detail="no Open Targets data")
        )

    # --- Clinical validation (Open Targets known_drug_score) -----------------
    if disease_assoc is not None:
        kd = disease_assoc.known_drug_score or 0.0
        if kd >= 0.7:
            rating, note = "strong", "approved / clinical-stage therapeutics exist"
        elif kd >= 0.3:
            rating, note = "moderate", "clinical-stage drug evidence"
        elif kd > 0:
            rating, note = "weak", "limited known-drug evidence"
        else:
            rating, note = "none", "no known-drug evidence"
        axes.append(
            EvidenceAxis(
                name="Clinical validation",
                rating=rating,
                detail=f"known-drug {kd:.2f} — {note}",
                value=kd,
            )
        )
    else:
        axes.append(
            EvidenceAxis(name="Clinical validation", rating="n/a", detail="no Open Targets data")
        )

    # --- Genetic evidence (Open Targets genetic_association is authoritative) --
    # OT's genetic_association_score already aggregates GWAS + colocalization +
    # curated genetics, so it ANCHORS this axis. GWAS Catalog hits are supporting
    # detail: a strong trait-relevant set can lift a moderate OT signal to strong,
    # but GWAS never independently inflates or downgrades an OT-backed verdict — so
    # brittle GWAS trait-label matching can't distort the rating. GWAS drives the
    # rating only when there is no OT data at all. The trait-relevance filter stays
    # a correctness check on which hits count. A high OT genetic score is the
    # monogenic/Mendelian signature, for which an empty GWAS Catalog is expected
    # rather than weak evidence.
    ot_genetic = (disease_assoc.genetic_association_score or 0.0) if disease_assoc else 0.0
    trait_relevant_hits = 0
    gwas_note: str | None = None
    if gwas_ev is not None:
        trait_relevant_hits = _trait_relevant_gwas_count(gwas_ev, indication)
        if trait_relevant_hits > 0:
            gwas_note = f"{trait_relevant_hits} trait-relevant GWAS hit(s)"
        elif gwas_ev.total_associations > 0:
            gwas_note = f"{gwas_ev.total_associations} GWAS hits but trait labels off-indication"
        else:
            gwas_note = "no GWAS hits"

    detail_parts: list[str] = []
    if disease_assoc is not None:
        detail_parts.append(f"OT genetic {ot_genetic:.2f}")
    if gwas_note:
        detail_parts.append(gwas_note)

    if gwas_ev is None and disease_assoc is None:
        axes.append(
            EvidenceAxis(
                name="Genetic evidence", rating="n/a", detail="no GWAS or Open Targets data"
            )
        )
    else:
        if disease_assoc is not None:
            # OT present → OT anchors the rating; trait-relevant GWAS corroborates upward only.
            if ot_genetic >= OT_GENETIC_MONOGENIC_THRESHOLD:
                rating = "strong"
                if trait_relevant_hits == 0:
                    detail_parts.append(
                        "monogenic/Mendelian signature — GWAS Catalog absent by design"
                    )
            elif ot_genetic >= 0.3:
                rating = "strong" if trait_relevant_hits >= 3 else "moderate"
            elif ot_genetic > 0:
                rating = "moderate" if trait_relevant_hits >= 1 else "weak"
            elif trait_relevant_hits >= 3:
                rating = "moderate"
            elif trait_relevant_hits >= 1:
                rating = "weak"
            else:
                rating = "none"
        else:
            # No OT data → GWAS Catalog is the only genetic signal and drives the rating.
            if trait_relevant_hits >= 3:
                rating = "strong"
            elif trait_relevant_hits >= 1:
                rating = "moderate"
            elif gwas_ev is not None and gwas_ev.total_associations > 0:
                rating = "weak"
            else:
                rating = "none"
        axes.append(
            EvidenceAxis(
                name="Genetic evidence",
                rating=rating,
                detail="; ".join(detail_parts) or "no genetic signal",
                value=ot_genetic or None,
            )
        )

    # --- Cancer dependency (DepMap CRISPR; OT somatic as proxy) --------------
    if cancer_dep is None:
        axes.append(EvidenceAxis(name="Cancer dependency", rating="n/a", detail="no DepMap data"))
    else:
        frac = cancer_dep.fraction_dependent_lines
        pct = int(frac * 100)
        is_real = "DepMap Chronos" in cancer_dep.data_source
        if cancer_dep.pan_essential:
            axes.append(
                EvidenceAxis(
                    name="Cancer dependency",
                    rating="weak",
                    detail=f"pan-essential ({pct}% of lines) — narrow therapeutic window",
                    value=frac,
                )
            )
        else:
            lineage_match = bool(indication) and any(
                indication.lower() in lin.lower() or lin.lower() in indication.lower()
                for lin in (cancer_dep.top_dependent_lineages or [])
            )
            rating = (
                "strong"
                if frac >= 0.5
                else "moderate"
                if frac >= 0.2
                else "weak"
                if frac > 0
                else "none"
            )
            # A low pan-cancer fraction concentrated in the queried lineage is
            # still strong evidence for *that* indication (qualitative, not a multiplier).
            if lineage_match and rating in ("moderate", "weak"):
                rating = "strong"
            # Proxy data (OT somatic, not real CRISPR) is lower confidence — cap at moderate.
            if not is_real and rating == "strong":
                rating = "moderate"
            parts = [f"{pct}% of lines dependent"]
            if not is_real:
                parts.append("OT somatic proxy")
            top = ", ".join((cancer_dep.top_dependent_lineages or [])[:3])
            if top:
                parts.append(f"top: {top}")
            if lineage_match:
                parts.append("matches indication")
            axes.append(
                EvidenceAxis(
                    name="Cancer dependency",
                    rating=rating,
                    detail="; ".join(parts),
                    value=frac,
                )
            )

    # --- Chemical tractability (ChEMBL potency, then PubChem breadth) --------
    # Functional/cell-based potency is rated a full tier above binding-only at
    # the same pChEMBL: binding ≠ cellular activity ≠ druggability.
    if chembl_compounds is not None and chembl_compounds.best_pchembl is not None:
        n = chembl_compounds.total_active_compounds
        bp_func = chembl_compounds.best_pchembl_functional
        bp_overall = chembl_compounds.best_pchembl
        if bp_func is not None:
            rating = "strong" if bp_func >= 8 else "moderate" if bp_func >= 6 else "weak"
            grade = (
                "clinical-grade"
                if bp_func >= 9
                else "lead-quality"
                if bp_func >= 7
                else "hit-quality"
            )
            detail = f"best functional pChEMBL {bp_func:.1f} ({grade}), {n} ChEMBL actives"
            value = bp_func
        else:
            rating = "moderate" if bp_overall >= 8 else "weak"
            detail = (
                f"best pChEMBL {bp_overall:.1f} (binding-only, no functional confirmation), "
                f"{n} ChEMBL actives"
            )
            value = bp_overall
        axes.append(
            EvidenceAxis(name="Chemical tractability", rating=rating, detail=detail, value=value)
        )
    elif compounds is not None and compounds.total_active_compounds >= PUBCHEM_MIN_COMPOUNDS:
        axes.append(
            EvidenceAxis(
                name="Chemical tractability",
                rating="weak",
                detail=f"{compounds.total_active_compounds} PubChem actives (count only, no potency data)",
                value=float(compounds.total_active_compounds),
            )
        )
    elif compounds is not None:
        axes.append(
            EvidenceAxis(
                name="Chemical tractability",
                rating="none",
                detail=(
                    f"{compounds.total_active_compounds} PubChem actives — below tractability "
                    f"threshold ({PUBCHEM_MIN_COMPOUNDS})"
                ),
                value=float(compounds.total_active_compounds),
            )
        )
    else:
        axes.append(
            EvidenceAxis(
                name="Chemical tractability", rating="n/a", detail="no ChEMBL or PubChem data"
            )
        )

    # --- Annotation quality (UniProt) ----------------------------------------
    if protein is None:
        axes.append(EvidenceAxis(name="Annotation quality", rating="n/a", detail="no UniProt data"))
    else:
        nvar = len(protein.known_variants)
        var_str = f"{nvar} known variant{'s' if nvar != 1 else ''}"
        if protein.reviewed:
            axes.append(
                EvidenceAxis(
                    name="Annotation quality",
                    rating="strong",
                    detail=f"SwissProt-reviewed; {var_str}",
                )
            )
        else:
            axes.append(
                EvidenceAxis(
                    name="Annotation quality",
                    rating="moderate",
                    detail=f"unreviewed (TrEMBL); {var_str}",
                )
            )

    # --- Tissue specificity (HPA) — only when extended mode fetched it -------
    # Restricted expression is a therapeutic-window advantage. Axis is added only
    # when HPA was queried; absent in standard mode rather than rated n/a.
    if protein_atlas is not None:
        cat = (
            protein_atlas.expression.rna_tissue_specificity_category
            if protein_atlas.expression
            else None
        )
        rating = _EXPRESSION_RATING.get(cat or "", "n/a")
        axes.append(
            EvidenceAxis(
                name="Tissue specificity",
                rating=rating,
                detail=f"HPA: {cat}" if cat else "no HPA specificity category",
            )
        )

    return axes


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
            match_note = (
                " — concentrated in the queried indication's lineage" if lineage_match else ""
            )
            parts.append(
                f"DepMap CRISPR data show dependency in {pct}% of cancer lines"
                + (f", highest in {top}" if top else "")
                + match_note
                + "."
            )

    gwas_relevant = _trait_relevant_gwas_count(gwas_ev, indication)
    if gwas_relevant > 0:
        p_str = f"{gwas_ev.strongest_p_value:.2e}" if gwas_ev.strongest_p_value else "N/A"
        parts.append(
            f"GWAS Catalog links {gwas_relevant} variant(s) near {symbol} "
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
        if gwas_relevant >= 5 and known_drug_score < 0.1 and literature_score < 0.15:
            ot_score = disease_assoc.overall_score if disease_assoc else 0.0
            parts.append(
                f"Caution: {gwas_relevant} GWAS hits near {symbol} for {indication} "
                f"with no known-drug evidence and sparse literature support "
                f"(OT drug score={known_drug_score:.2f}, literature={literature_score:.2f}, "
                f"overall={ot_score:.2f}). This pattern is consistent with LD with a nearby causal "
                f"variant, comorbidity confounding, or indirect drug effects — not direct target biology. "
                f"Functional validation is recommended before treating GWAS signal as target evidence."
            )
    elif gwas_ev and gwas_ev.total_associations > 0:
        parts.append(
            f"GWAS Catalog returns {gwas_ev.total_associations} hit(s) near {symbol}, but their "
            f"trait labels do not match '{indication}' — not counted as indication genetic evidence."
        )

    if chembl_compounds and chembl_compounds.best_pchembl is not None:
        bp = chembl_compounds.best_pchembl
        ic50_nm = 10 ** (9 - bp)
        # Use models._format_ic50_nm so sub-nM potency renders as "50 pM" not "0.0 nM" (Bug Q).
        from genesis_bio_mcp.models import _format_ic50_nm

        ic50_str = _format_ic50_nm(ic50_nm)
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
