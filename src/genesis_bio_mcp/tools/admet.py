"""Local ADMET / developability screening — pure RDKit, no network, no API key.

A *keyless* first tier for the ADMET capability gap the eval harness flagged. Everything
here is a rule-based heuristic computed from RDKit descriptors — **screening alerts, not
validated predictions**:

- aqueous solubility via the Delaney **ESOL** linear model,
- GI-absorption / BBB heuristics from TPSA + logP,
- a coarse **hERG** flag (basic amine + high lipophilicity — the classic pharmacophore),
- a coarse **CYP3A4** metabolic-liability flag (lipophilic, higher-MW scaffold),
- **Brenk + PAINS** structural alerts,
- Ertl **SAscore** synthetic accessibility (ships with RDKit Contrib).

Quantitative ML endpoints (hERG/CYP IC50, clearance, t½) are an explicit future optional
tier (ADMET-AI / TDC) — kept off the keyless serving path by design.
"""

from __future__ import annotations

import logging
import os
import sys

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, RDConfig
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

from genesis_bio_mcp.models import ADMETProfile
from genesis_bio_mcp.tools.cheminformatics import _standardized_parent

logger = logging.getLogger(__name__)

# SAscore ships in RDKit's Contrib tree, not the importable package — add it to the path once.
_SA_PATH = os.path.join(RDConfig.RDContribDir, "SA_Score")
if _SA_PATH not in sys.path:
    sys.path.append(_SA_PATH)
try:
    import sascorer  # type: ignore[import-not-found]
except Exception:  # pragma: no cover - defensive: SAscore missing in some RDKit builds
    sascorer = None
    logger.warning(
        "RDKit SA_Score contrib not available; synthesizability will be reported Unknown"
    )

# Brenk (unwanted functionality) + PAINS (assay interference) structural-alert catalogs.
_ALERT_PARAMS = FilterCatalogParams()
_ALERT_PARAMS.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
_ALERT_PARAMS.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
_ALERT_CATALOG = FilterCatalog(_ALERT_PARAMS)

# A protonatable aliphatic amine (basic nitrogen): trivalent N that is not an amide,
# imine, charged, or aromatic-attached (anilines are only weakly basic). The basic-centre
# + lipophilicity combination is the canonical hERG-binding pharmacophore.
_BASIC_AMINE = Chem.MolFromSmarts("[NX3;H0,H1,H2;!$(NC=[O,S,N]);!$(N=*);!$([N+]);!$(N-a)]")


def _esol_logs(mol, mw: float, logp: float) -> float:
    """Delaney (2004) ESOL estimate of aqueous log S (mol/L)."""
    rb = Descriptors.NumRotatableBonds(mol)
    heavy = mol.GetNumHeavyAtoms() or 1
    aromatic_proportion = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic()) / heavy
    logs = 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rb - 0.74 * aromatic_proportion
    return round(logs, 2)


def _solubility_class(logs: float) -> str:
    if logs >= -2:
        return "highly soluble"
    if logs >= -4:
        return "soluble"
    if logs >= -6:
        return "poorly soluble"
    return "insoluble"


def _herg(mol, logp: float) -> tuple[str, str]:
    has_basic_amine = _BASIC_AMINE is not None and mol.HasSubstructMatch(_BASIC_AMINE)
    if has_basic_amine and logp >= 3.5:
        return (
            "Elevated",
            f"basic amine + high lipophilicity (logP {logp:.1f}) — classic hERG pharmacophore",
        )
    if logp >= 4.5:
        return "Moderate", f"high lipophilicity (logP {logp:.1f})"
    if has_basic_amine and logp >= 2.5:
        return "Moderate", f"basic amine with moderate lipophilicity (logP {logp:.1f})"
    return "Low", "no basic-amine / high-logP hERG pharmacophore"


def _cyp3a4(mw: float, logp: float, aromatic_rings: int) -> tuple[str, str]:
    if logp >= 3.0 and mw >= 350 and aromatic_rings >= 2:
        return (
            "Elevated",
            "lipophilic, higher-MW polyaromatic scaffold — typical CYP3A4 substrate space",
        )
    if logp >= 4.0:
        return (
            "Elevated",
            f"high lipophilicity (logP {logp:.1f}) favors CYP3A4 oxidative metabolism",
        )
    return "Low", "not in the lipophilic high-MW CYP3A4 substrate space"


def assess_admet(smiles: str) -> ADMETProfile | None:
    """Compute the rule-based ADMET / developability panel for a SMILES, or None if unparseable."""
    if not smiles or not smiles.strip():
        return None
    cleaned = smiles.strip()
    mol = Chem.MolFromSmiles(cleaned)
    if mol is None:
        logger.debug("RDKit could not parse SMILES for ADMET: %s", cleaned)
        return None

    parent = _standardized_parent(mol)  # salt-strip + neutralize so logP/MW reflect the drug
    mw = Descriptors.MolWt(parent)
    logp = Crippen.MolLogP(parent)
    tpsa = Descriptors.TPSA(parent)
    aromatic_rings = Descriptors.NumAromaticRings(parent)

    logs = _esol_logs(parent, mw, logp)
    # Veber-style oral-absorption heuristic: low polar surface area + drug-like lipophilicity.
    gi_absorption = "High" if (tpsa <= 140 and -1.0 <= logp <= 5.0) else "Low"
    bbb_permeant = tpsa <= 90 and mw <= 450 and logp >= 1.0
    herg_risk, herg_reason = _herg(parent, logp)
    cyp_liability, cyp_reason = _cyp3a4(mw, logp, aromatic_rings)

    alerts = sorted({m.GetDescription() for m in _ALERT_CATALOG.GetMatches(parent)})

    sa_score = None
    synthesizability = "Unknown"
    if sascorer is not None:
        try:
            sa_score = round(sascorer.calculateScore(parent), 2)
            synthesizability = (
                "Easy" if sa_score <= 3.0 else "Moderate" if sa_score <= 6.0 else "Difficult"
            )
        except Exception as exc:  # pragma: no cover - defensive
            logger.warning("SAscore failed for %s: %s", cleaned, exc)

    return ADMETProfile(
        input_smiles=cleaned,
        canonical_smiles=Chem.MolToSmiles(parent),
        molecular_weight=round(mw, 2),
        logp=round(logp, 2),
        tpsa=round(tpsa, 2),
        esol_logs=logs,
        solubility_class=_solubility_class(logs),
        gi_absorption=gi_absorption,
        bbb_permeant=bbb_permeant,
        herg_risk=herg_risk,
        herg_reason=herg_reason,
        cyp3a4_liability=cyp_liability,
        cyp3a4_reason=cyp_reason,
        structural_alerts=alerts,
        sa_score=sa_score,
        synthesizability=synthesizability,
    )
