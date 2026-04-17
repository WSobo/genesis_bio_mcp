"""Pure-Python protein biochemistry + liability-motif scanner.

Used by :mod:`genesis_bio_mcp.server.get_protein_sequence` and (planned)
the antibody CDR developability extension of :mod:`sabdab`.

No Biopython dependency — the underlying formulas are small and
well-documented. Each constant block below cites its source so a reviewer
can sanity-check the numbers against the published reference.

References
----------
* Average amino-acid residue masses: ExPASy ProtParam (Gasteiger et al.,
  2005) average isotopic masses.
* pKa values: Bjellqvist et al., Electrophoresis 14:1023, 1993 — the set
  used by ExPASy ProtParam. We use a simplified uniform-N-term variant
  (position-specific N-term pKas are a small ~0.05 pI refinement that
  would double the table size for little gain on the reference proteins).
* Hydropathy scale: Kyte & Doolittle, J. Mol. Biol. 157:105, 1982.
* Molar extinction coefficient at 280 nm: Edelhoch, Biochemistry 6:1948,
  1967 (revised by Pace et al., Protein Sci. 4:2411, 1995) —
  W=5500, Y=1490, cystine (disulfide, S-S)=125 M⁻¹cm⁻¹ in H₂O.
"""

from __future__ import annotations

import logging
import re
from typing import Literal

from pydantic import BaseModel, Field

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Amino-acid constants — each value documented with source + rationale.
# ---------------------------------------------------------------------------

# Standard 20-letter amino-acid alphabet (one-letter codes).
_AMINO_ACIDS: frozenset[str] = frozenset("ACDEFGHIKLMNPQRSTVWY")

# Average residue masses (Da), per ExPASy ProtParam. Add 18.01528 Da (one
# water) per complete chain to account for the free N-/C-terminus.
# Source: https://web.expasy.org/protparam/protparam-doc.html
_RESIDUE_MASSES_DA: dict[str, float] = {
    "A": 71.0788,
    "R": 156.1875,
    "N": 114.1038,
    "D": 115.0886,
    "C": 103.1388,
    "E": 129.1155,
    "Q": 128.1307,
    "G": 57.0519,
    "H": 137.1411,
    "I": 113.1594,
    "L": 113.1594,
    "K": 128.1741,
    "M": 131.1926,
    "F": 147.1766,
    "P": 97.1167,
    "S": 87.0782,
    "T": 101.1051,
    "W": 186.2132,
    "Y": 163.1760,
    "V": 99.1326,
}
_WATER_MASS_DA: float = 18.01528

# Bjellqvist pKa values — the set used by ExPASy ProtParam.
# Source: Bjellqvist et al., Electrophoresis 14:1023, 1993; refined in
# Bjellqvist et al., Electrophoresis 15:529, 1994. Our reference-protein
# tests confirm this set matches ProtParam pI values to within ±0.1 pH.
_PK_N_TERM: float = 7.5
_PK_C_TERM: float = 3.55
_PK_SIDECHAIN: dict[str, float] = {
    # Positively charged at neutral pH → protonated form is the charged form.
    "K": 10.0,
    "R": 12.0,
    "H": 5.98,
    # Negatively charged at neutral pH → deprotonated form is charged.
    "D": 4.05,
    "E": 4.45,
    "C": 9.0,  # thiol; ionization produces thiolate (-1)
    "Y": 10.0,  # phenol; ionization produces phenolate (-1)
}
_POSITIVE_SIDECHAINS: frozenset[str] = frozenset("KRH")
_NEGATIVE_SIDECHAINS: frozenset[str] = frozenset("DECY")

# Kyte & Doolittle hydropathy index. GRAVY = mean over residues.
# Tuning note: values are immutable — this is the reference scale.
_KD_HYDROPATHY: dict[str, float] = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}

# Extinction coefficients at 280 nm (M⁻¹ cm⁻¹) — Edelhoch / Pace revision.
# Cystine (oxidized disulfide) contributes; free cysteine does not absorb
# appreciably at 280 nm.
_EXT_W: int = 5500
_EXT_Y: int = 1490
_EXT_CYSTINE: int = 125  # per S-S bond (2 cysteines)

_AROMATIC: frozenset[str] = frozenset("FWY")

# Liability-motif regexes. Positions reported are 1-indexed, pointing to the
# residue that *starts* the motif. Documented with the mechanistic reason
# the motif is flagged.
#
# Tuning note: widening these patterns (e.g. including DH for isomerization
# or QG for deamidation) trades recall for false-positive rate.  The set
# below follows the consensus hotspot table in Xu et al., MAbs 11(2):239,
# 2019, which is the standard reference for therapeutic-protein developability.
_LIABILITY_PATTERNS: list[tuple[str, str, str]] = [
    # (regex, motif_type, mechanism_description)
    (r"N[GS]", "deamidation", "Asn→Asp/isoAsp under neutral/basic pH"),
    (r"D[GS]", "isomerization", "Asp→isoAsp via succinimide intermediate"),
    (r"N[^P][ST]", "n_glycosylation", "N-X-S/T sequon (X≠P); glycan attachment site"),
]


# ---------------------------------------------------------------------------
# Output model
# ---------------------------------------------------------------------------


class LiabilityHit(BaseModel):
    """One hit from the liability-motif scanner."""

    motif_type: Literal[
        "deamidation",
        "isomerization",
        "n_glycosylation",
        "oxidation_methionine",
        "oxidation_tryptophan",
        "free_cysteine",
        "cysteine_position",
    ] = Field(description="Liability category — see biochem module docstring for mechanisms")
    position: int = Field(ge=1, description="1-indexed residue position where the motif starts")
    residues: str = Field(description="Matched residues at the hit site")
    context: str = Field(
        default="",
        description="Up to ±3 residues around the hit for quick inspection",
    )


class BiochemFeatures(BaseModel):
    """Computed biochemical features of a protein sequence."""

    length: int = Field(ge=0, description="Residue count")
    molecular_weight_Da: float = Field(description="Average molecular weight in Daltons")
    theoretical_pI: float = Field(
        description="Theoretical isoelectric point (pH at which net charge = 0)"
    )
    gravy: float = Field(
        description="Grand average of hydropathy (Kyte-Doolittle); positive = hydrophobic"
    )
    net_charge_pH74: float = Field(description="Net charge at physiological pH 7.4")
    aromatic_fraction: float = Field(
        ge=0.0, le=1.0, description="Fraction of residues that are F/W/Y"
    )
    cysteine_count: int = Field(ge=0, description="Total cysteine residues")
    extinction_coefficient_280nm: int = Field(
        description="ε₂₈₀ assuming all Cys form disulfides (oxidized form, M⁻¹ cm⁻¹)"
    )
    extinction_coefficient_280nm_reduced: int = Field(
        description="ε₂₈₀ assuming all Cys are free thiols (reduced form, M⁻¹ cm⁻¹)"
    )


# ---------------------------------------------------------------------------
# Sequence sanitization
# ---------------------------------------------------------------------------


def _clean(seq: str) -> str:
    """Uppercase, strip whitespace, and drop characters outside the 20-AA alphabet.

    Non-standard codes (X, B, Z, U, *, -) are silently dropped with a debug
    log. This matches ProtParam's behavior — unknown residues should not
    poison biochem calculations.
    """
    raw = seq.upper().strip()
    cleaned = "".join(c for c in raw if c in _AMINO_ACIDS)
    if len(cleaned) != len(raw.replace(" ", "").replace("\n", "").replace("\t", "")):
        logger.debug(
            "biochem: dropped %d non-standard residues from %d-aa input",
            len(raw) - len(cleaned),
            len(raw),
        )
    return cleaned


# ---------------------------------------------------------------------------
# Core biochem functions
# ---------------------------------------------------------------------------


def molecular_weight(seq: str) -> float:
    """Average molecular weight in Daltons."""
    s = _clean(seq)
    if not s:
        return 0.0
    return sum(_RESIDUE_MASSES_DA[aa] for aa in s) + _WATER_MASS_DA


def _net_charge_at_ph(seq: str, ph: float) -> float:
    """Net charge of the chain at a given pH via Henderson-Hasselbalch."""
    if not seq:
        return 0.0
    # N-terminus: protonated form carries +1
    pos = 1.0 / (1.0 + 10.0 ** (ph - _PK_N_TERM))
    # C-terminus: deprotonated form carries -1
    neg = 1.0 / (1.0 + 10.0 ** (_PK_C_TERM - ph))
    # Side chains
    for aa in seq:
        pka = _PK_SIDECHAIN.get(aa)
        if pka is None:
            continue
        if aa in _POSITIVE_SIDECHAINS:
            pos += 1.0 / (1.0 + 10.0 ** (ph - pka))
        else:
            neg += 1.0 / (1.0 + 10.0 ** (pka - ph))
    return pos - neg


def net_charge(seq: str, ph: float = 7.4) -> float:
    """Net charge at given pH (default 7.4, physiological)."""
    return _net_charge_at_ph(_clean(seq), ph)


def theoretical_pi(seq: str, tol: float = 1e-4, max_iter: int = 100) -> float:
    """Theoretical isoelectric point via bisection on net charge over pH ∈ [0, 14].

    Convergence: iterate until ``|net_charge| < tol`` or ``max_iter`` reached.
    The function is monotonically decreasing in pH on the well-posed domain,
    so bisection is guaranteed to converge for any non-empty chain.

    Raises:
        RuntimeError: if ``max_iter`` is exhausted without convergence.
            This indicates a bug (numerical instability or a pathological
            sequence) rather than a legitimate input — log and bubble up.

    Tuning note: ``tol=1e-4`` matches ProtParam output precision (pI reported
    to 2 decimal places). Loosening tol to 1e-2 speeds up by ~6 iterations
    at the cost of a potential ±0.01 pH drift.
    """
    s = _clean(seq)
    if not s:
        return 7.0
    low, high = 0.0, 14.0
    for _ in range(max_iter):
        mid = (low + high) / 2.0
        q = _net_charge_at_ph(s, mid)
        if abs(q) < tol:
            return mid
        if q > 0:
            low = mid
        else:
            high = mid
    raise RuntimeError(
        f"theoretical_pi did not converge in {max_iter} iterations; "
        f"last net_charge={q:.6f} at pH={mid:.4f}"
    )


def gravy(seq: str) -> float:
    """Kyte-Doolittle grand average of hydropathy."""
    s = _clean(seq)
    if not s:
        return 0.0
    return sum(_KD_HYDROPATHY[aa] for aa in s) / len(s)


def aromatic_fraction(seq: str) -> float:
    """Fraction of residues that are F, W, or Y."""
    s = _clean(seq)
    if not s:
        return 0.0
    return sum(1 for aa in s if aa in _AROMATIC) / len(s)


def cysteine_positions(seq: str) -> list[int]:
    """1-indexed positions of every Cys residue."""
    s = _clean(seq)
    return [i + 1 for i, aa in enumerate(s) if aa == "C"]


def extinction_coefficient_280nm(seq: str, *, reduced: bool = True) -> int:
    """Molar extinction coefficient at 280 nm (M⁻¹ cm⁻¹), Edelhoch formula.

    Args:
        reduced: If True (default), assume all cysteines are free thiols
            (their 280 nm absorption is negligible). If False, assume all
            cysteines are paired in disulfides, contributing ε_cystine per
            S-S bond. Real proteins are mixtures — report both and let
            the caller choose.
    """
    s = _clean(seq)
    if not s:
        return 0
    n_w = s.count("W")
    n_y = s.count("Y")
    eps = n_w * _EXT_W + n_y * _EXT_Y
    if not reduced:
        n_c = s.count("C")
        eps += (n_c // 2) * _EXT_CYSTINE
    return eps


def compute_features(seq: str) -> BiochemFeatures:
    """Compute the full BiochemFeatures bundle for a sequence in one call."""
    s = _clean(seq)
    return BiochemFeatures(
        length=len(s),
        molecular_weight_Da=round(molecular_weight(s), 2),
        theoretical_pI=round(theoretical_pi(s), 2),
        gravy=round(gravy(s), 3),
        net_charge_pH74=round(net_charge(s, 7.4), 2),
        aromatic_fraction=round(aromatic_fraction(s), 3),
        cysteine_count=s.count("C"),
        extinction_coefficient_280nm=extinction_coefficient_280nm(s, reduced=False),
        extinction_coefficient_280nm_reduced=extinction_coefficient_280nm(s, reduced=True),
    )


# ---------------------------------------------------------------------------
# Liability scanner
# ---------------------------------------------------------------------------


def _context(seq: str, position_1idx: int, span: int = 3) -> str:
    """Return up to ±span residues around a 1-indexed position for display."""
    i = position_1idx - 1
    lo = max(0, i - span)
    hi = min(len(seq), i + span + 1)
    prefix = "…" if lo > 0 else ""
    suffix = "…" if hi < len(seq) else ""
    return f"{prefix}{seq[lo:hi]}{suffix}"


def scan_liabilities(
    seq: str,
    disulfide_annotated_positions: set[int] | None = None,
) -> list[LiabilityHit]:
    """Scan a sequence for common protein-engineering liability motifs.

    Liabilities scanned:
      * **Deamidation hotspots** — NG, NS (Asn→Asp/isoAsp over shelf-life).
      * **Isomerization hotspots** — DG, DS (Asp→isoAsp via succinimide).
      * **N-glycosylation sequons** — N-X-S/T where X≠P (glycan attachment).
      * **Oxidation-prone residues** — every M and W (reported as data, not
        flagged as an immediate risk — these are common and context-dependent).
      * **Cysteine positions** — every Cys residue. If
        ``disulfide_annotated_positions`` is provided (from UniProt DISULFID
        features), Cys residues NOT in that set are flagged as
        ``free_cysteine`` (a real liability for aggregation / scrambling).
        Without annotation, Cys positions are reported as
        ``cysteine_position`` data points without a liability label —
        parity-based guessing (odd total count) is not reliable.

    Args:
        seq: Protein sequence (1-letter codes).
        disulfide_annotated_positions: Optional set of 1-indexed Cys positions
            known to be disulfide-bonded per UniProt DISULFID features.

    Returns:
        List of :class:`LiabilityHit`, sorted by position ascending.
    """
    s = _clean(seq)
    hits: list[LiabilityHit] = []
    if not s:
        return hits

    for pattern, motif_type, _reason in _LIABILITY_PATTERNS:
        for m in re.finditer(pattern, s):
            pos = m.start() + 1
            hits.append(
                LiabilityHit(
                    motif_type=motif_type,  # type: ignore[arg-type]
                    position=pos,
                    residues=m.group(0),
                    context=_context(s, pos),
                )
            )

    for i, aa in enumerate(s):
        if aa == "M":
            pos = i + 1
            hits.append(
                LiabilityHit(
                    motif_type="oxidation_methionine",
                    position=pos,
                    residues="M",
                    context=_context(s, pos),
                )
            )
        elif aa == "W":
            pos = i + 1
            hits.append(
                LiabilityHit(
                    motif_type="oxidation_tryptophan",
                    position=pos,
                    residues="W",
                    context=_context(s, pos),
                )
            )
        elif aa == "C":
            pos = i + 1
            if disulfide_annotated_positions is not None:
                if pos not in disulfide_annotated_positions:
                    hits.append(
                        LiabilityHit(
                            motif_type="free_cysteine",
                            position=pos,
                            residues="C",
                            context=_context(s, pos),
                        )
                    )
            else:
                hits.append(
                    LiabilityHit(
                        motif_type="cysteine_position",
                        position=pos,
                        residues="C",
                        context=_context(s, pos),
                    )
                )

    hits.sort(key=lambda h: h.position)
    return hits
