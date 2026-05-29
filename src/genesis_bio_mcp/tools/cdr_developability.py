"""Antibody CDR developability assessment.

Pure-Python analysis layered on :mod:`genesis_bio_mcp.tools.biochem`. Given the
six CDR sequences of an antibody (auto-numbered from VH/VL via AbNum, or supplied
directly), it computes per-CDR biochemical features and scans for the
developability liabilities that matter most in the binding loops:

* Sequence-motif liabilities (deamidation NG/NS, isomerization DG/DS,
  N-glycosylation sequons, Met/Trp oxidation, cysteines) — via
  :func:`biochem.scan_liabilities`.
* Property risks — high hydrophobicity (aggregation / poly-specificity) and an
  unusually long CDR-H3 (developability / PK risk).

CDRs are where these liabilities are most consequential: motif chemistry there
can erode binding over shelf-life, and CDR-H3 hydrophobicity/length drive
aggregation and clearance. Thresholds below cite the consensus therapeutic-
antibody developability literature (Raybould et al., PNAS 116:4025, 2019 — the
"Therapeutic Antibody Profiler"; Xu et al., MAbs 11(2):239, 2019).
"""

from __future__ import annotations

from pydantic import BaseModel, Field

from genesis_bio_mcp.tools.biochem import LiabilityHit, compute_features, scan_liabilities

# Canonical CDR field order + display labels (Chothia field names from sabdab).
_CDR_ORDER: list[str] = ["vh_cdr1", "vh_cdr2", "vh_cdr3", "vl_cdr1", "vl_cdr2", "vl_cdr3"]
_CDR_LABELS: dict[str, str] = {
    "vh_cdr1": "VH CDR1",
    "vh_cdr2": "VH CDR2",
    "vh_cdr3": "VH CDR3",
    "vl_cdr1": "VL CDR1",
    "vl_cdr2": "VL CDR2",
    "vl_cdr3": "VL CDR3",
}

# Property thresholds. Tuning notes:
#  * GRAVY > 0 already means net-hydrophobic; >0.5 in a CDR is a meaningful
#    aggregation / poly-specificity flag per TAP-style profiling.
#  * CDR-H3 ≥ 18 aa sits in the long tail of the human repertoire and
#    correlates with developability / PK liabilities.
_GRAVY_HYDROPHOBIC = 0.5
_LONG_CDRH3 = 18

# Motif types that warrant an explicit risk flag (vs. informational only).
_FLAGGED_MOTIFS: dict[str, str] = {
    "deamidation": "deamidation motif",
    "isomerization": "isomerization motif",
    "n_glycosylation": "N-glycosylation sequon",
    "oxidation_tryptophan": "oxidation-prone Trp",
}
_CYS_MOTIFS = frozenset({"free_cysteine", "cysteine_position"})


class CDRLiabilitySummary(BaseModel):
    """Per-CDR biochemical + liability summary."""

    cdr: str = Field(description="CDR field key, e.g. 'vh_cdr3'")
    label: str = Field(description="Human-readable CDR label, e.g. 'VH CDR3'")
    sequence: str = Field(description="CDR amino-acid sequence")
    length: int = Field(description="CDR length in residues")
    net_charge_pH74: float = Field(description="Net charge of the CDR at pH 7.4")
    gravy: float = Field(description="Kyte-Doolittle GRAVY (positive = hydrophobic)")
    liabilities: list[LiabilityHit] = Field(
        default_factory=list,
        description="Liability-motif hits within this CDR (positions are CDR-local)",
    )


class CDRDevelopabilityReport(BaseModel):
    """Developability assessment across an antibody's CDRs."""

    numbering_source: str = Field(
        description="How CDRs were obtained, e.g. 'AbNum (Chothia)' or 'user-provided'"
    )
    cdrs: list[CDRLiabilitySummary] = Field(
        default_factory=list, description="Per-CDR summaries in canonical order"
    )
    total_liabilities: int = Field(
        default=0, description="Total liability-motif hits across all CDRs"
    )
    risk_flags: list[str] = Field(
        default_factory=list,
        description="Human-readable high-priority developability callouts",
    )

    def to_markdown(self) -> str:
        lines = [
            "## Antibody CDR Developability",
            f"_CDR source: {self.numbering_source}_",
        ]
        if not self.cdrs:
            lines.append(
                "\n_No CDR sequences could be determined. Provide VH/VL sequences "
                "(for AbNum numbering) or explicit CDR sequences._"
            )
            return "\n".join(lines)

        lines += [
            "",
            "| CDR | Sequence | Len | Charge | GRAVY | Liabilities |",
            "|---|---|---|---|---|---|",
        ]
        for c in self.cdrs:
            lines.append(
                f"| {c.label} | `{c.sequence}` | {c.length} | {c.net_charge_pH74:+.1f} "
                f"| {c.gravy:+.2f} | {len(c.liabilities)} |"
            )

        if self.risk_flags:
            lines += ["", "### Developability flags", ""]
            lines += [f"- {flag}" for flag in self.risk_flags]
        else:
            lines += ["", "_No high-priority developability liabilities flagged in the CDRs._"]

        lines += [
            "",
            f"**Total liability-motif hits across CDRs:** {self.total_liabilities}",
            "",
            "_Motif chemistry (deamidation/isomerization/glycosylation/oxidation) in CDRs can "
            "erode binding over shelf-life; CDR-H3 hydrophobicity and length drive aggregation "
            "and clearance. Flags are heuristics — confirm against structure and assay data._",
        ]
        return "\n".join(lines)


def assess_cdr_developability(
    cdr_map: dict[str, str], numbering_source: str
) -> CDRDevelopabilityReport:
    """Build a :class:`CDRDevelopabilityReport` from a map of CDR field → sequence.

    Args:
        cdr_map: Mapping of CDR field keys (``'vh_cdr1'`` … ``'vl_cdr3'``, plus any
            extra keys) to amino-acid sequences. Empty / missing CDRs are skipped.
        numbering_source: Provenance string shown in the report.

    Returns:
        A populated :class:`CDRDevelopabilityReport` (empty ``cdrs`` if nothing usable).
    """
    summaries: list[CDRLiabilitySummary] = []
    risk_flags: list[str] = []
    total = 0

    ordered_keys = _CDR_ORDER + [k for k in cdr_map if k not in _CDR_ORDER]
    for key in ordered_keys:
        raw = (cdr_map.get(key) or "").strip()
        if not raw:
            continue
        feats = compute_features(raw)
        if feats.length == 0:
            continue
        liabs = scan_liabilities(raw)
        total += len(liabs)
        label = _CDR_LABELS.get(key, key)
        summaries.append(
            CDRLiabilitySummary(
                cdr=key,
                label=label,
                sequence=raw.upper(),
                length=feats.length,
                net_charge_pH74=feats.net_charge_pH74,
                gravy=feats.gravy,
                liabilities=liabs,
            )
        )

        for h in liabs:
            if h.motif_type in _FLAGGED_MOTIFS:
                risk_flags.append(
                    f"{label}: {_FLAGGED_MOTIFS[h.motif_type]} ({h.residues}) at CDR position {h.position}"
                )
            elif h.motif_type in _CYS_MOTIFS:
                risk_flags.append(
                    f"{label}: cysteine at CDR position {h.position} — potential unpaired thiol "
                    "(disulfide scrambling / aggregation risk)"
                )
        if feats.gravy > _GRAVY_HYDROPHOBIC:
            risk_flags.append(
                f"{label}: high hydrophobicity (GRAVY={feats.gravy:+.2f}) — aggregation / "
                "poly-specificity risk"
            )
        if key == "vh_cdr3" and feats.length >= _LONG_CDRH3:
            risk_flags.append(
                f"{label}: long CDR-H3 ({feats.length} aa) — may affect developability / PK"
            )

    return CDRDevelopabilityReport(
        numbering_source=numbering_source,
        cdrs=summaries,
        total_liabilities=total,
        risk_flags=risk_flags,
    )
