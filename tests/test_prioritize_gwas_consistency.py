"""Regression: the prioritize_target prose summary must not contradict the
per-axis genetic rating about GWAS support.

`_build_evidence_profile` applies a trait-relevance gate — off-indication GWAS
hits are labelled "off-indication" and not counted as genetic evidence. The
prose summary used to ignore that gate and assert "GWAS Catalog links N variants
... to '{indication}'-related traits" for ANY raw hit count, so a target could
read "GWAS off-indication" in its axis table and "GWAS supports" in its summary
at the same time. Both now route through `_trait_relevant_gwas_count`.
"""

from __future__ import annotations

from genesis_bio_mcp.models import GwasEvidence, GwasHit
from genesis_bio_mcp.tools.target_prioritization import (
    _build_summary,
    _trait_relevant_gwas_count,
)


def _hit(trait: str) -> GwasHit:
    return GwasHit(
        study_accession="GCST000001",
        trait=trait,
        mapped_gene="GENE",
        risk_allele="rs1-A",
        p_value=1e-9,
    )


def _gwas(*, n: int, trait: str, trait_query: str) -> GwasEvidence:
    return GwasEvidence(
        gene_symbol="GENE",
        trait_query=trait_query,
        total_associations=n,
        associations=[_hit(trait) for _ in range(min(n, 3))],
        strongest_p_value=1e-9,
    )


def test_trait_relevant_count_gate() -> None:
    assert _trait_relevant_gwas_count(None, "melanoma") == 0
    # Off-indication labels -> 0 (token-set Jaccard below threshold).
    assert (
        _trait_relevant_gwas_count(
            _gwas(n=8, trait="Immunoglobulin levels", trait_query="melanoma"), "melanoma"
        )
        == 0
    )
    # Fallback (non-exact trait query) -> 0 regardless of labels.
    assert (
        _trait_relevant_gwas_count(
            _gwas(n=8, trait="melanoma", trait_query="melanoma (no exact-trait match)"), "melanoma"
        )
        == 0
    )
    # On-indication labels -> full count.
    assert (
        _trait_relevant_gwas_count(_gwas(n=4, trait="melanoma", trait_query="melanoma"), "melanoma")
        == 4
    )


def test_summary_does_not_claim_off_indication_gwas_support() -> None:
    off = _gwas(n=8, trait="Immunoglobulin levels", trait_query="melanoma")
    md = _build_summary("GENE", "melanoma", None, None, off, None)
    # The false positive claim must be gone...
    assert "to 'melanoma'-related traits" not in md
    assert "links 8" not in md
    # ...and the honest off-indication framing present, matching the axis.
    assert "do not match 'melanoma'" in md
    assert "not counted as indication genetic evidence" in md


def test_summary_reports_trait_relevant_gwas_support() -> None:
    on = _gwas(n=4, trait="melanoma", trait_query="melanoma")
    md = _build_summary("GENE", "melanoma", None, None, on, None)
    assert "GWAS Catalog links 4 variant(s) near GENE" in md
    assert "to 'melanoma'-related traits" in md
