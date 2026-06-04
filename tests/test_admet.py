"""Tests for the local ADMET / developability panel (tools/admet.py).

Pure RDKit compute — no network, no key. The heuristics are coarse by design;
these tests pin behavior on well-known molecules and the structural contract,
not exact physicochemical values.
"""

from __future__ import annotations

from genesis_bio_mcp.models import ADMETProfile
from genesis_bio_mcp.tools.admet import assess_admet

# A kinase inhibitor with a basic (piperazine) amine + lipophilic, polyaromatic scaffold —
# clinically a known hERG binder and CYP3A4 substrate.
IMATINIB = "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"


def test_invalid_smiles_returns_none():
    assert assess_admet("not_a_smiles") is None
    assert assess_admet("") is None
    assert assess_admet("   ") is None


def test_imatinib_flags_herg_and_cyp_liabilities():
    p = assess_admet(IMATINIB)
    assert isinstance(p, ADMETProfile)
    assert p.herg_risk == "Elevated"  # basic amine + high logP
    assert p.cyp3a4_liability == "Elevated"  # lipophilic, higher-MW polyaromatic
    assert p.molecular_weight > 400
    assert p.sa_score is not None and p.synthesizability in {"Easy", "Moderate", "Difficult"}


def test_aspirin_is_low_liability_and_soluble():
    p = assess_admet(ASPIRIN)
    assert p.herg_risk == "Low"
    assert p.cyp3a4_liability == "Low"
    # Small, polar acid → ESOL should not call it insoluble.
    assert p.esol_logs > -4
    assert p.solubility_class in {"highly soluble", "soluble"}


def test_salt_is_stripped_before_scoring():
    """A salt form should score like its parent (sodium acetate ~ acetic acid)."""
    parent = assess_admet("CC(=O)O")
    salt = assess_admet("CC(=O)[O-].[Na+]")
    assert parent is not None and salt is not None
    assert salt.canonical_smiles == parent.canonical_smiles
    assert abs(salt.molecular_weight - parent.molecular_weight) < 0.01


def test_markdown_surfaces_herg_cyp_and_synthesizability():
    md = assess_admet(IMATINIB).to_markdown()
    assert "hERG" in md
    assert "CYP3A4" in md
    assert "SAscore" in md
    # Honesty caveat must be present.
    assert "not validated predictions" in md.lower()


def test_structural_alerts_is_a_list():
    p = assess_admet(IMATINIB)
    assert isinstance(p.structural_alerts, list)
