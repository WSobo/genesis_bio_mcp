"""Unit tests for the evidence-profile engine (`_build_evidence_profile`).

These replace the old `_compute_score` / `ScoreBreakdown` tests. Instead of pinning
a composite number for named anecdote targets, each test checks that one axis is
*classified* correctly from its source metric (strong / moderate / weak / none / n-a).
Boundaries are conventional cutoffs, applied independently per axis — nothing is
summed or weighted into a single score.
"""

from genesis_bio_mcp.models import (
    CancerDependency,
    ChEMBLCompounds,
    ComparisonReport,
    Compounds,
    EvidenceAxis,
    GeneResolution,
    GwasEvidence,
    GwasHit,
    HPAExpression,
    ProteinAtlasReport,
    ProteinInfo,
    TargetComparisonRow,
    TargetDiseaseAssociation,
    TargetPrioritizationReport,
)
from genesis_bio_mcp.tools.target_prioritization import _build_evidence_profile

# ---------------------------------------------------------------------------
# Builders for synthetic inputs
# ---------------------------------------------------------------------------


def _assoc(**kw) -> TargetDiseaseAssociation:
    base = dict(
        gene_symbol="X",
        disease_name="d",
        disease_efo_id="EFO_1",
        ensembl_id="ENSG1",
        overall_score=0.5,
        evidence_count=10,
    )
    base.update(kw)
    return TargetDiseaseAssociation(**base)


def _dep(**kw) -> CancerDependency:
    base = dict(
        gene_symbol="X",
        mean_ceres_score=-0.5,
        fraction_dependent_lines=0.6,
        pan_essential=False,
        top_dependent_lineages=["Skin"],
        cell_lines=[],
        data_source="DepMap Chronos Combined (300/1100 cell lines dependent)",
    )
    base.update(kw)
    return CancerDependency(**base)


def _gwas(n: int, trait: str = "melanoma", *, fallback: bool = False) -> GwasEvidence:
    tq = (
        f"{trait} (no exact-trait match — top gene-level associations shown)" if fallback else trait
    )
    hits = [
        GwasHit(
            study_accession=f"GCST{i}",
            trait=trait,
            mapped_gene="X",
            risk_allele=f"rs{i}-A",
            p_value=1e-12,
        )
        for i in range(n)
    ]
    return GwasEvidence(
        gene_symbol="X",
        trait_query=tq,
        total_associations=n,
        associations=hits,
        strongest_p_value=1e-12 if n else None,
    )


def _chembl(**kw) -> ChEMBLCompounds:
    base = dict(gene_symbol="X", total_active_compounds=50, compounds=[])
    base.update(kw)
    return ChEMBLCompounds(**base)


def _profile(**kw) -> dict[str, EvidenceAxis]:
    """Run the engine and return {axis_name: axis}."""
    return {a.name: a for a in _build_evidence_profile(**kw)}


def _bare(**kw):
    """_build_evidence_profile with all sources None unless overridden."""
    args = dict(
        disease_assoc=None,
        cancer_dep=None,
        gwas_ev=None,
        compounds=None,
        protein=None,
        chembl_compounds=None,
        indication="melanoma",
    )
    args.update(kw)
    return _profile(**args)


# ---------------------------------------------------------------------------
# Disease association axis
# ---------------------------------------------------------------------------


def test_disease_association_boundaries():
    assert _bare(disease_assoc=_assoc(overall_score=0.5))["Disease association"].rating == "strong"
    assert (
        _bare(disease_assoc=_assoc(overall_score=0.2))["Disease association"].rating == "moderate"
    )
    assert _bare(disease_assoc=_assoc(overall_score=0.05))["Disease association"].rating == "weak"
    assert _bare(disease_assoc=_assoc(overall_score=0.0))["Disease association"].rating == "none"
    assert _bare(disease_assoc=None)["Disease association"].rating == "n/a"


# ---------------------------------------------------------------------------
# Clinical validation axis
# ---------------------------------------------------------------------------


def test_clinical_validation_boundaries():
    assert (
        _bare(disease_assoc=_assoc(known_drug_score=0.7))["Clinical validation"].rating == "strong"
    )
    assert (
        _bare(disease_assoc=_assoc(known_drug_score=0.3))["Clinical validation"].rating
        == "moderate"
    )
    assert _bare(disease_assoc=_assoc(known_drug_score=0.1))["Clinical validation"].rating == "weak"
    assert _bare(disease_assoc=_assoc(known_drug_score=0.0))["Clinical validation"].rating == "none"
    assert _bare(disease_assoc=None)["Clinical validation"].rating == "n/a"


# ---------------------------------------------------------------------------
# Genetic evidence axis
# ---------------------------------------------------------------------------


def test_genetic_evidence_monogenic_is_strong_without_gwas():
    """OT genetic > 0.7 (Mendelian signature) → strong even with no GWAS hits."""
    axis = _bare(disease_assoc=_assoc(genetic_association_score=0.99), gwas_ev=None)[
        "Genetic evidence"
    ]
    assert axis.rating == "strong"
    assert "monogenic" in axis.detail.lower()


def test_genetic_evidence_polygenic_without_gwas_is_moderate():
    """OT genetic 0.55 (below monogenic threshold), no GWAS → moderate, not strong."""
    axis = _bare(disease_assoc=_assoc(genetic_association_score=0.55), gwas_ev=None)[
        "Genetic evidence"
    ]
    assert axis.rating == "moderate"


def test_genetic_evidence_three_trait_relevant_hits_is_strong():
    axis = _bare(gwas_ev=_gwas(3, "melanoma"), indication="melanoma")["Genetic evidence"]
    assert axis.rating == "strong"


def test_genetic_evidence_two_trait_relevant_hits_is_moderate():
    axis = _bare(gwas_ev=_gwas(2, "melanoma"), indication="melanoma")["Genetic evidence"]
    assert axis.rating == "moderate"


def test_genetic_evidence_fallback_hits_not_credited():
    """Fallback (no exact-trait match) hits are off-trait → weak, never strong."""
    axis = _bare(gwas_ev=_gwas(5, "melanoma", fallback=True), indication="melanoma")[
        "Genetic evidence"
    ]
    assert axis.rating == "weak"


def test_genetic_evidence_off_trait_hits_not_credited():
    """Direct-match hits whose trait labels don't overlap the indication → weak."""
    off = GwasEvidence(
        gene_symbol="TP53",
        trait_query="li-fraumeni syndrome",
        total_associations=3,
        associations=[
            GwasHit(
                study_accession="GCST7",
                trait="sex hormone-binding globulin",
                mapped_gene="TP53",
                risk_allele="rs9-A",
                p_value=4e-276,
            )
        ],
        strongest_p_value=4e-276,
    )
    axis = _bare(gwas_ev=off, indication="li-fraumeni syndrome")["Genetic evidence"]
    assert axis.rating == "weak"


def test_genetic_evidence_na_when_no_sources():
    assert _bare(gwas_ev=None, disease_assoc=None)["Genetic evidence"].rating == "n/a"


# ---------------------------------------------------------------------------
# Cancer dependency axis
# ---------------------------------------------------------------------------


def test_cancer_dependency_boundaries():
    assert (
        _bare(cancer_dep=_dep(fraction_dependent_lines=0.6))["Cancer dependency"].rating == "strong"
    )
    assert (
        _bare(cancer_dep=_dep(fraction_dependent_lines=0.3, top_dependent_lineages=["Blood"]))[
            "Cancer dependency"
        ].rating
        == "moderate"
    )
    assert _bare(cancer_dep=None)["Cancer dependency"].rating == "n/a"


def test_cancer_dependency_pan_essential_is_weak():
    axis = _bare(cancer_dep=_dep(pan_essential=True, fraction_dependent_lines=0.95))[
        "Cancer dependency"
    ]
    assert axis.rating == "weak"
    assert "pan-essential" in axis.detail.lower()


def test_cancer_dependency_proxy_caps_at_moderate():
    """OT-somatic proxy (no real CRISPR) cannot be rated strong even at high fraction."""
    proxy = _dep(
        fraction_dependent_lines=0.8,
        data_source="Open Targets somatic mutation (approximated proxy)",
        top_dependent_lineages=["Blood"],
    )
    axis = _bare(cancer_dep=proxy, indication="leukemia")["Cancer dependency"]
    assert axis.rating == "moderate"
    assert "proxy" in axis.detail.lower()


def test_cancer_dependency_lineage_match_lifts_low_fraction():
    """A low pan-cancer fraction concentrated in the queried lineage → strong."""
    dep = _dep(fraction_dependent_lines=0.1, top_dependent_lineages=["melanoma"])
    axis = _bare(cancer_dep=dep, indication="melanoma")["Cancer dependency"]
    assert axis.rating == "strong"
    assert "matches indication" in axis.detail


# ---------------------------------------------------------------------------
# Chemical tractability axis
# ---------------------------------------------------------------------------


def test_chemical_tractability_functional_vs_binding():
    # Functional clinical-grade → strong
    assert (
        _bare(chembl_compounds=_chembl(best_pchembl=9.2, best_pchembl_functional=9.2))[
            "Chemical tractability"
        ].rating
        == "strong"
    )
    # Functional lead-quality (6–8) → moderate
    assert (
        _bare(chembl_compounds=_chembl(best_pchembl=7.5, best_pchembl_functional=7.5))[
            "Chemical tractability"
        ].rating
        == "moderate"
    )
    # Binding-only clinical-grade → moderate (one tier below functional)
    assert (
        _bare(
            chembl_compounds=_chembl(
                best_pchembl=9.2, best_pchembl_functional=None, best_pchembl_binding=9.2
            )
        )["Chemical tractability"].rating
        == "moderate"
    )
    # Binding-only lead-quality (<8) → weak
    assert (
        _bare(
            chembl_compounds=_chembl(
                best_pchembl=7.5, best_pchembl_functional=None, best_pchembl_binding=7.5
            )
        )["Chemical tractability"].rating
        == "weak"
    )


def test_chemical_tractability_pubchem_fallback():
    # PubChem count above threshold, no ChEMBL → weak
    assert (
        _bare(compounds=Compounds(gene_symbol="X", total_active_compounds=50, compounds=[]))[
            "Chemical tractability"
        ].rating
        == "weak"
    )
    # Below tractability threshold → none
    assert (
        _bare(compounds=Compounds(gene_symbol="X", total_active_compounds=2, compounds=[]))[
            "Chemical tractability"
        ].rating
        == "none"
    )
    # No chemistry data at all → n/a
    assert _bare(compounds=None, chembl_compounds=None)["Chemical tractability"].rating == "n/a"


# ---------------------------------------------------------------------------
# Annotation quality axis
# ---------------------------------------------------------------------------


def _protein(reviewed: bool) -> ProteinInfo:
    return ProteinInfo(
        uniprot_accession="P1",
        gene_symbol="X",
        protein_name="p",
        organism="Homo sapiens",
        function_summary="",
        subcellular_locations=[],
        pathways=[],
        disease_associations=[],
        pdb_structures=[],
        known_variants=[],
        reviewed=reviewed,
    )


def test_annotation_quality():
    assert _bare(protein=_protein(True))["Annotation quality"].rating == "strong"
    assert _bare(protein=_protein(False))["Annotation quality"].rating == "moderate"
    assert _bare(protein=None)["Annotation quality"].rating == "n/a"


# ---------------------------------------------------------------------------
# Tissue specificity axis (extended mode only)
# ---------------------------------------------------------------------------


def _atlas(category) -> ProteinAtlasReport:
    return ProteinAtlasReport(
        gene_symbol="X",
        expression=HPAExpression(gene_symbol="X", rna_tissue_specificity_category=category),
        pathology=[],
    )


def test_tissue_specificity_rating_map():
    assert _bare(protein_atlas=_atlas("Tissue enriched"))["Tissue specificity"].rating == "strong"
    assert _bare(protein_atlas=_atlas("Tissue enhanced"))["Tissue specificity"].rating == "strong"
    assert _bare(protein_atlas=_atlas("Group enriched"))["Tissue specificity"].rating == "moderate"
    assert (
        _bare(protein_atlas=_atlas("Low tissue specificity"))["Tissue specificity"].rating == "weak"
    )
    assert _bare(protein_atlas=_atlas("Not detected"))["Tissue specificity"].rating == "none"
    assert _bare(protein_atlas=_atlas("Bogus category"))["Tissue specificity"].rating == "n/a"


def test_tissue_specificity_absent_in_standard_mode():
    """Without HPA (standard mode), the tissue axis is omitted entirely — not rated n/a."""
    assert "Tissue specificity" not in _bare(protein_atlas=None)


# ---------------------------------------------------------------------------
# Report rendering
# ---------------------------------------------------------------------------


def test_report_renders_evidence_profile_not_score():
    report = TargetPrioritizationReport(
        gene_symbol="BRAF",
        indication="melanoma",
        resolution=GeneResolution(hgnc_symbol="BRAF", source="input"),
        disease_association=_assoc(overall_score=0.89, known_drug_score=0.88),
        evidence_profile=_build_evidence_profile(
            _assoc(overall_score=0.89, known_drug_score=0.88),
            None,
            None,
            None,
            None,
            indication="melanoma",
        ),
        evidence_summary="BRAF strong in melanoma.",
    )
    md = report.to_markdown()
    assert "## Evidence Profile" in md
    assert "Strong" in md
    # No composite score / tier anywhere
    assert "/10" not in md
    assert "Priority Score" not in md
    assert "Disease association" in md


# ---------------------------------------------------------------------------
# Comparison matrix ordering
# ---------------------------------------------------------------------------


def _row(gene: str, ratings: dict[str, str]) -> TargetComparisonRow:
    return TargetComparisonRow(
        gene_symbol=gene,
        evidence_profile=[
            EvidenceAxis(name=name, rating=r, detail="d") for name, r in ratings.items()
        ],
        evidence_summary=f"{gene} summary",
    )


def test_comparison_matrix_orders_by_dominance_and_disclaims_score():
    strong_target = _row(
        "AAA",
        {
            "Disease association": "strong",
            "Cancer dependency": "strong",
            "Chemical tractability": "strong",
        },
    )
    weak_target = _row(
        "ZZZ",
        {
            "Disease association": "weak",
            "Cancer dependency": "none",
            "Chemical tractability": "weak",
        },
    )
    # Pass weak first to prove ordering isn't input order.
    md = ComparisonReport(indication="melanoma", rows=[weak_target, strong_target]).to_markdown()

    # Strong target ranks first.
    assert md.index("AAA") < md.index("ZZZ")
    # Honest framing + no composite score.
    assert "not" in md.lower() and "validated composite score" in md
    assert "/10" not in md
    # Matrix uses the short axis labels and glyph legend.
    assert "Legend:" in md


def test_comparison_matrix_breaks_ties_alphabetically():
    a = _row("BBB", {"Disease association": "strong"})
    b = _row("AAA", {"Disease association": "strong"})
    md = ComparisonReport(indication="x", rows=[a, b]).to_markdown()
    assert md.index("AAA") < md.index("BBB")
