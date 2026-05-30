"""Local cheminformatics — RDKit-computed molecular properties.

Pure, offline, deterministic compute (no network, no rate limits). Turns a SMILES
string into physicochemical descriptors and drug-likeness assessments used across
small-molecule discovery: MW, logP, TPSA, H-bond donors/acceptors, rotatable bonds,
QED, Lipinski/Veber rule checks, a PAINS assay-interference flag, and the
Bemis-Murcko scaffold.

This is the first ``Methods4Insight``-style tool: it *computes* rather than
retrieves, complementing the database lookups (PubChem/ChEMBL) with structure-level
analysis.
"""

from __future__ import annotations

import logging

from rdkit import Chem
from rdkit.Chem import QED, Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

from genesis_bio_mcp.models import MolecularProperties, StandardizedStructure

logger = logging.getLogger(__name__)

# PAINS (Pan-Assay INterference compoundS) substructure catalog — built once.
# A match flags a molecule prone to false-positive assay readouts.
_PAINS_PARAMS = FilterCatalogParams()
_PAINS_PARAMS.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
_PAINS_CATALOG = FilterCatalog(_PAINS_PARAMS)


def compute_molecular_properties(smiles: str) -> MolecularProperties | None:
    """Compute physicochemical + drug-likeness properties from a SMILES string.

    Returns ``None`` if RDKit cannot parse the SMILES (invalid structure).
    """
    if not smiles or not smiles.strip():
        return None
    cleaned = smiles.strip()
    mol = Chem.MolFromSmiles(cleaned)
    if mol is None:
        logger.debug("RDKit could not parse SMILES: %s", cleaned)
        return None

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = Descriptors.TPSA(mol)
    rot = Descriptors.NumRotatableBonds(mol)

    # Lipinski Rule of Five: each breached threshold is one violation; ≤1 is
    # conventionally still considered drug-like.
    violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])

    # Bemis-Murcko scaffold; empty string for acyclic molecules → None.
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol) or None

    return MolecularProperties(
        input_smiles=cleaned,
        canonical_smiles=Chem.MolToSmiles(mol),
        molecular_formula=rdMolDescriptors.CalcMolFormula(mol),
        molecular_weight=round(mw, 2),
        logp=round(logp, 2),
        tpsa=round(tpsa, 2),
        h_bond_donors=hbd,
        h_bond_acceptors=hba,
        rotatable_bonds=rot,
        aromatic_rings=Lipinski.NumAromaticRings(mol),
        heavy_atom_count=mol.GetNumHeavyAtoms(),
        fraction_csp3=round(Descriptors.FractionCSP3(mol), 3),
        qed=round(QED.qed(mol), 3),
        lipinski_violations=violations,
        passes_lipinski=violations <= 1,
        passes_veber=(rot <= 10 and tpsa <= 140),
        pains_alert=_PAINS_CATALOG.HasMatch(mol),
        murcko_scaffold=scaffold,
    )


def standardize_structure(smiles: str) -> StandardizedStructure | None:
    """Normalize a molecule: strip salts/solvents, neutralize charges, canonical tautomer.

    Produces a registration-ready canonical form plus InChI / InChIKey — the
    data-readiness primitive for deduplicating and joining compounds across sources.
    Returns ``None`` if RDKit cannot parse the SMILES.
    """
    if not smiles or not smiles.strip():
        return None
    cleaned = smiles.strip()
    mol = Chem.MolFromSmiles(cleaned)
    if mol is None:
        logger.debug("RDKit could not parse SMILES for standardization: %s", cleaned)
        return None

    try:
        sanitized = rdMolStandardize.Cleanup(mol)
        # Largest organic fragment — strips salt counter-ions and solvents.
        parent = rdMolStandardize.FragmentParent(sanitized)
        # Neutralize charges where chemically sensible.
        neutral = rdMolStandardize.Uncharger().uncharge(parent)
        # Canonical tautomer for consistent registration.
        canonical = rdMolStandardize.TautomerEnumerator().Canonicalize(neutral)
    except Exception as exc:  # RDKit can raise on pathological inputs
        logger.warning("RDKit standardization failed for %s: %s", cleaned, exc)
        return None

    std_smiles = Chem.MolToSmiles(canonical)
    try:
        inchi = Chem.MolToInchi(canonical) or None
        inchikey = Chem.MolToInchiKey(canonical) or None
    except Exception:
        inchi, inchikey = None, None

    return StandardizedStructure(
        input_smiles=cleaned,
        standardized_smiles=std_smiles,
        inchi=inchi,
        inchikey=inchikey,
        molecular_formula=rdMolDescriptors.CalcMolFormula(canonical),
        salt_or_solvent_removed=canonical.GetNumAtoms() < sanitized.GetNumAtoms(),
        changed=std_smiles != Chem.MolToSmiles(mol),
    )
