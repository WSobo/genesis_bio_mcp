"""Unit tests for the RDKit-based cheminformatics tools."""

from genesis_bio_mcp.models import BatchMolecularProperties
from genesis_bio_mcp.tools.cheminformatics import (
    compute_molecular_properties,
    standardize_structure,
)


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


def test_standardize_strips_salt_and_neutralizes():
    # Sodium acetate → acetic acid (Na+ stripped, carboxylate neutralized).
    s = standardize_structure("CC(=O)[O-].[Na+]")
    assert s is not None
    assert s.standardized_smiles == "CC(=O)O"
    assert s.molecular_formula == "C2H4O2"
    assert s.salt_or_solvent_removed is True
    assert s.changed is True
    assert s.inchikey == "QTBSBXVTEAMEQO-UHFFFAOYSA-N"  # acetic acid
    md = s.to_markdown()
    assert "Standardized Structure" in md
    assert "salt/solvent stripped" in md


def test_standardize_noop_on_clean_molecule():
    # Aspirin is already standard — no salt, no change.
    s = standardize_structure("CC(=O)Oc1ccccc1C(=O)O")
    assert s is not None
    assert s.standardized_smiles == "CC(=O)Oc1ccccc1C(=O)O"
    assert s.salt_or_solvent_removed is False
    assert s.changed is False


def test_standardize_invalid_smiles_returns_none():
    assert standardize_structure("not_a_smiles!!!") is None
    assert standardize_structure("") is None


def test_batch_molecular_properties_markdown_summary():
    # Two valid molecules + one parse failure → one summary row each, failures listed.
    valid = [
        compute_molecular_properties("CCO"),
        compute_molecular_properties("CC(=O)Oc1ccccc1C(=O)O"),
    ]
    batch = BatchMolecularProperties(
        results=[p for p in valid if p is not None],
        failed=["not_a_smiles!!!"],
        total_requested=3,
    )
    md = batch.to_markdown()
    assert "3 molecules" in md
    assert "2 computed, 1 failed" in md
    # One table row per computed molecule
    assert "C9H8O4" in md  # aspirin formula
    assert "C2H6O" in md  # ethanol formula
    assert "Failed to parse (1)" in md
    assert "`not_a_smiles!!!`" in md


def test_batch_molecular_properties_all_failed_has_no_table():
    batch = BatchMolecularProperties(results=[], failed=["x", "y"], total_requested=2)
    md = batch.to_markdown()
    assert "0 computed, 2 failed" in md
    assert "| # | Input |" not in md  # no property table when nothing parsed
    assert "Failed to parse (2)" in md
