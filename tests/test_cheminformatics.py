"""Unit tests for the RDKit-based cheminformatics tool."""

from genesis_bio_mcp.tools.cheminformatics import compute_molecular_properties


def test_compute_properties_aspirin():
    p = compute_molecular_properties("CC(=O)Oc1ccccc1C(=O)O")
    assert p is not None
    assert p.molecular_formula == "C9H8O4"
    assert round(p.molecular_weight) == 180
    assert p.h_bond_donors == 1
    assert p.h_bond_acceptors == 3
    assert p.aromatic_rings == 1
    assert p.rotatable_bonds == 2
    assert p.passes_lipinski is True
    assert p.passes_veber is True
    assert p.pains_alert is False
    assert p.murcko_scaffold == "c1ccccc1"
    md = p.to_markdown()
    assert "C9H8O4" in md
    assert "Lipinski" in md
    assert "Bemis-Murcko" in md


def test_compute_properties_canonicalizes_input():
    # A non-canonical SMILES for aspirin should canonicalize identically.
    p = compute_molecular_properties("O=C(O)c1ccccc1OC(C)=O")
    assert p is not None
    assert p.canonical_smiles == "CC(=O)Oc1ccccc1C(=O)O"


def test_compute_properties_invalid_smiles_returns_none():
    assert compute_molecular_properties("not_a_smiles!!!") is None
    assert compute_molecular_properties("") is None
    assert compute_molecular_properties("   ") is None


def test_compute_properties_acyclic_has_no_scaffold():
    # Ethanol is acyclic → Murcko scaffold is empty → None.
    p = compute_molecular_properties("CCO")
    assert p is not None
    assert p.murcko_scaffold is None
    assert p.aromatic_rings == 0


def test_compute_properties_flags_non_druglike_molecule():
    # A C40 alkane breaks two Lipinski rules (MW > 500 and logP > 5) → fails Ro5.
    p = compute_molecular_properties("C" * 40)
    assert p is not None
    assert p.lipinski_violations >= 2
    assert p.passes_lipinski is False
