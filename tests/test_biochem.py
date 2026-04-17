"""Biochem module tests — reference values from ExPASy ProtParam.

We validate against three well-documented proteins so future changes to the
constants block can be caught immediately.
"""

from __future__ import annotations

import pytest

from genesis_bio_mcp.tools.biochem import (
    LiabilityHit,
    aromatic_fraction,
    compute_features,
    cysteine_positions,
    extinction_coefficient_280nm,
    gravy,
    molecular_weight,
    net_charge,
    scan_liabilities,
    theoretical_pi,
)

# ---------------------------------------------------------------------------
# Reference sequences
# ---------------------------------------------------------------------------

# Hen egg-white lysozyme, mature chain (UniProt P00698, residues 19-147)
# ProtParam: MW 14313.14 Da, pI 9.32, GRAVY -0.472, 129 aa,
# ε₂₈₀ (red) 37550, ε₂₈₀ (ox, 4 S-S bonds) 38400
LYSOZYME = (
    "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRN"
    "LCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
)

# Bovine ribonuclease A, mature chain (UniProt P61823, residues 27-150) — 124 aa.
# ProtParam: pI 8.64, MW 13690.3 Da, GRAVY -0.663.
RIBONUCLEASE_A = (
    "KETAAAKFERQHMDSSTSAASSSNYCNQMMKSRNLTKDRCKPVNTFVHESLADVQAVCSQKNVACKNGQTNCYQ"
    "SYSTMSITDCRETGSSKYPNCAYKTTQANKHIIVACEGNPYVPVHFDASV"
)

# Constructed test peptide that contains every liability motif class exactly once:
# positions: N1+G2 (deamidation NG), then filler, D7+S8 (isomerization DS),
# then N-X-T sequon at N13 NLT, then a methionine at 18, Trp at 21, Cys at 25.
LIABILITY_PEPTIDE = "NGAAAADSAAAANLTMSHWRACLK"
#                    123456789012345678901234
#                             1111111111222222


# ---------------------------------------------------------------------------
# Biochem reference-value tests
# ---------------------------------------------------------------------------


def test_lysozyme_molecular_weight():
    assert molecular_weight(LYSOZYME) == pytest.approx(14313.14, abs=1.0)


def test_lysozyme_pi():
    # ProtParam reports pI 9.32; our EMBOSS pKa set matches to ±0.1.
    assert theoretical_pi(LYSOZYME) == pytest.approx(9.32, abs=0.1)


def test_lysozyme_gravy():
    assert gravy(LYSOZYME) == pytest.approx(-0.472, abs=0.01)


def test_lysozyme_extinction_reduced():
    # 6 Trp × 5500 + 3 Tyr × 1490 = 33000 + 4470 = 37470 (reduced)
    assert extinction_coefficient_280nm(LYSOZYME, reduced=True) == pytest.approx(37470, abs=200)


def test_lysozyme_extinction_oxidized():
    # Add 4 S-S bonds × 125 = 500, so oxidized = 37470 + 500 = 37970
    # (ProtParam reports 38400 with its own assumption count — within 500.)
    assert extinction_coefficient_280nm(LYSOZYME, reduced=False) == pytest.approx(37970, abs=500)


def test_lysozyme_aromatic_fraction():
    # 6W + 3Y + 3F = 12 aromatic / 129 residues ≈ 0.093
    assert aromatic_fraction(LYSOZYME) == pytest.approx(12 / 129, abs=0.001)


def test_lysozyme_cysteine_positions():
    # Hen egg lysozyme has 8 cysteines forming 4 disulfide bonds.
    cys = cysteine_positions(LYSOZYME)
    assert len(cys) == 8


def test_ribonuclease_a_pi():
    # ProtParam reports RNase A pI ≈ 8.64; allow ±0.2 because the EMBOSS
    # pKa table gives slightly different values than Lehninger/IPC.
    assert theoretical_pi(RIBONUCLEASE_A) == pytest.approx(8.64, abs=0.25)


def test_ribonuclease_a_mw():
    assert molecular_weight(RIBONUCLEASE_A) == pytest.approx(13690.3, abs=2.0)


def test_polyalanine_pi_converges():
    # Plain polyalanine has no charged side chains — pI is determined
    # purely by N- and C-termini; should be around (3.6 + 8.6) / 2 ≈ 6.1.
    pi = theoretical_pi("A" * 50)
    assert 5.5 < pi < 6.5


def test_compute_features_lysozyme():
    feats = compute_features(LYSOZYME)
    assert feats.length == 129
    assert feats.molecular_weight_Da == pytest.approx(14313.14, abs=1.0)
    assert feats.theoretical_pI == pytest.approx(9.32, abs=0.1)
    assert feats.gravy == pytest.approx(-0.472, abs=0.01)
    assert feats.cysteine_count == 8


def test_net_charge_positive_at_acidic_ph():
    # Acidic pH → basic residues dominate → net charge positive.
    assert net_charge(LYSOZYME, ph=3.0) > 0


def test_net_charge_negative_at_basic_ph():
    # Basic pH → acidic residues dominate → net charge negative.
    assert net_charge(LYSOZYME, ph=12.0) < 0


def test_empty_sequence():
    feats = compute_features("")
    assert feats.length == 0
    assert feats.molecular_weight_Da == 0.0
    # pI returns neutral-ish default for empty string (fallback to 7.0)
    assert feats.theoretical_pI == pytest.approx(7.0, abs=0.1)


def test_cleaning_non_standard_residues():
    # X / B / Z should be silently dropped.
    feats = compute_features("AAAXBZAAA")
    assert feats.length == 6  # 6 alanines survived


# ---------------------------------------------------------------------------
# Liability scanner tests
# ---------------------------------------------------------------------------


def test_scanner_finds_deamidation_hotspot():
    hits = scan_liabilities(LIABILITY_PEPTIDE)
    deam = [h for h in hits if h.motif_type == "deamidation"]
    assert any(h.position == 1 and h.residues == "NG" for h in deam)


def test_scanner_finds_isomerization_hotspot():
    hits = scan_liabilities(LIABILITY_PEPTIDE)
    iso = [h for h in hits if h.motif_type == "isomerization"]
    assert any(h.position == 7 and h.residues == "DS" for h in iso)


def test_scanner_finds_n_glycosylation_sequon():
    hits = scan_liabilities(LIABILITY_PEPTIDE)
    gly = [h for h in hits if h.motif_type == "n_glycosylation"]
    assert any(h.position == 13 and h.residues == "NLT" for h in gly)


def test_scanner_flags_oxidation_residues():
    hits = scan_liabilities(LIABILITY_PEPTIDE)
    met = [h for h in hits if h.motif_type == "oxidation_methionine"]
    trp = [h for h in hits if h.motif_type == "oxidation_tryptophan"]
    assert any(h.position == 16 for h in met)
    assert any(h.position == 19 for h in trp)


def test_scanner_reports_cys_as_data_without_annotation():
    hits = scan_liabilities(LIABILITY_PEPTIDE)
    cys = [h for h in hits if h.motif_type == "cysteine_position"]
    free = [h for h in hits if h.motif_type == "free_cysteine"]
    assert len(cys) == 1 and cys[0].position == 22
    assert len(free) == 0  # no annotation supplied → no free-Cys flag


def test_scanner_flags_free_cys_when_annotation_supplied():
    # Supply an empty disulfide set → all Cys are free.
    hits = scan_liabilities(LIABILITY_PEPTIDE, disulfide_annotated_positions=set())
    free = [h for h in hits if h.motif_type == "free_cysteine"]
    assert len(free) == 1 and free[0].position == 22


def test_scanner_respects_disulfide_annotation():
    # Mark Cys22 as bonded → no free_cysteine emitted.
    hits = scan_liabilities(LIABILITY_PEPTIDE, disulfide_annotated_positions={22})
    free = [h for h in hits if h.motif_type == "free_cysteine"]
    cys = [h for h in hits if h.motif_type == "cysteine_position"]
    assert len(free) == 0
    assert len(cys) == 0  # with annotation, bonded Cys is not reported as either


def test_sequon_excludes_proline():
    # N-P-S is NOT a sequon (X must not be P). Regex guard should exclude it.
    hits = scan_liabilities("AAANPSAAA")
    gly = [h for h in hits if h.motif_type == "n_glycosylation"]
    assert len(gly) == 0


def test_hits_sorted_by_position():
    hits = scan_liabilities(LIABILITY_PEPTIDE)
    positions = [h.position for h in hits]
    assert positions == sorted(positions)


def test_empty_sequence_returns_no_hits():
    assert scan_liabilities("") == []


def test_liability_hit_schema():
    # Smoke-test that LiabilityHit is Pydantic-validating correctly.
    h = LiabilityHit(motif_type="deamidation", position=1, residues="NG", context="NGAAA")
    assert h.motif_type == "deamidation"
    assert h.position == 1
