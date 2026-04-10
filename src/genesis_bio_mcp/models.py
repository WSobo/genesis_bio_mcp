"""Pydantic output models for all MCP tool responses."""

from __future__ import annotations

from typing import Optional
from pydantic import BaseModel, Field, model_validator


# ---------------------------------------------------------------------------
# Gene resolution
# ---------------------------------------------------------------------------


class GeneResolution(BaseModel):
    """Canonical gene identifiers resolved from a synonym or alias."""

    hgnc_symbol: str = Field(description="Approved HGNC gene symbol, e.g. 'BRAF'")
    hgnc_id: Optional[str] = Field(None, description="HGNC ID, e.g. 'HGNC:1097'")
    ncbi_gene_id: Optional[str] = Field(None, description="NCBI Entrez Gene ID")
    uniprot_accession: Optional[str] = Field(None, description="Primary UniProt accession")
    synonyms: list[str] = Field(default_factory=list, description="Known gene aliases")
    source: str = Field(description="Resolution source: 'uniprot' | 'ncbi' | 'input'")

    def to_markdown(self) -> str:
        parts = [f"**{self.hgnc_symbol}**"]
        if self.ncbi_gene_id:
            parts.append(f"NCBI Gene: {self.ncbi_gene_id}")
        if self.uniprot_accession:
            parts.append(f"UniProt: {self.uniprot_accession}")
        if self.hgnc_id:
            parts.append(f"HGNC: {self.hgnc_id}")
        line = " | ".join(parts)
        syns = f"\nAliases: {', '.join(self.synonyms[:8])}" if self.synonyms else ""
        return f"## Gene: {line}{syns}"


# ---------------------------------------------------------------------------
# UniProt / protein info
# ---------------------------------------------------------------------------


class KnownVariant(BaseModel):
    """A disease-linked natural variant annotated in UniProt."""

    position: Optional[str] = Field(None, description="Sequence position (e.g. '600')")
    original: Optional[str] = Field(None, description="Original amino acid (single letter)")
    variant: Optional[str] = Field(None, description="Variant amino acid (single letter)")
    disease: Optional[str] = Field(None, description="Associated disease name")
    clinical_significance: Optional[str] = Field(None, description="Pathogenicity annotation if available")


class ProteinInfo(BaseModel):
    """Protein-level annotation from UniProt Swiss-Prot."""

    uniprot_accession: str = Field(description="Primary UniProt accession, e.g. 'P15056'")
    gene_symbol: str = Field(description="HGNC gene symbol")
    protein_name: str = Field(description="Recommended protein name from UniProt")
    organism: str = Field(description="Source organism")
    function_summary: str = Field(description="Curated functional description (CC FUNCTION)")
    subcellular_locations: list[str] = Field(description="Annotated subcellular compartments")
    pathways: list[str] = Field(description="Reactome pathway names linked in UniProt")
    disease_associations: list[str] = Field(description="MIM disease entries linked in UniProt")
    pdb_structures: list[str] = Field(description="PDB accessions with experimental structures")
    known_variants: list[KnownVariant] = Field(
        default_factory=list, description="Curated disease-linked natural variants"
    )
    reviewed: bool = Field(description="True if Swiss-Prot (manually reviewed), False if TrEMBL")

    def to_markdown(self) -> str:
        reviewed_label = "Swiss-Prot (manually reviewed)" if self.reviewed else "TrEMBL (unreviewed)"
        lines = [
            f"## Protein: {self.protein_name} ({self.gene_symbol})",
            f"**UniProt:** {self.uniprot_accession} — {reviewed_label}",
            "",
            "### Function",
            self.function_summary,
        ]
        if self.subcellular_locations:
            lines += ["", f"**Subcellular location:** {', '.join(self.subcellular_locations)}"]
        if self.pathways:
            lines += ["", f"**Pathways ({len(self.pathways)}):** {', '.join(self.pathways[:5])}"]
            if len(self.pathways) > 5:
                lines[-1] += f" (+{len(self.pathways) - 5} more)"
        if self.disease_associations:
            lines += ["", f"**Disease associations:** {', '.join(self.disease_associations[:5])}"]
        if self.pdb_structures:
            lines += ["", f"**PDB structures ({len(self.pdb_structures)}):** {', '.join(self.pdb_structures[:8])}"]
        if self.known_variants:
            lines += ["", "### Disease-linked Variants"]
            for v in self.known_variants[:5]:
                pos = f"p.{v.original}{v.position}{v.variant}" if (v.original and v.position and v.variant) else v.position or "?"
                disease = v.disease or "unknown disease"
                sig = f" ({v.clinical_significance})" if v.clinical_significance else ""
                lines.append(f"- **{pos}** — {disease}{sig}")
            if len(self.known_variants) > 5:
                lines.append(f"- _...and {len(self.known_variants) - 5} more_")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Open Targets
# ---------------------------------------------------------------------------


class DiseaseLinkEvidence(BaseModel):
    """A single evidence type contributing to an Open Targets association."""

    evidence_type: str = Field(
        description="Evidence category: 'genetic_association' | 'somatic_mutation' | "
        "'known_drug' | 'literature' | 'animal_model' | 'rna_expression'"
    )
    score: float = Field(description="Evidence score 0–1 for this datatype")


class TargetDiseaseAssociation(BaseModel):
    """Open Targets evidence-based association between a gene target and a disease."""

    gene_symbol: str
    disease_name: str
    disease_efo_id: str = Field(description="EFO ontology ID resolved from disease name")
    ensembl_id: str = Field(description="Ensembl gene ID used for the query")
    overall_score: float = Field(ge=0.0, le=1.0, description="Aggregate association score 0–1; >0.5 is strong")
    genetic_association_score: Optional[float] = Field(None, description="Score from genetic evidence (GWAS, rare variants)")
    somatic_mutation_score: Optional[float] = Field(None, description="Score from somatic mutations in cancer")
    known_drug_score: Optional[float] = Field(None, description="Score from approved/clinical-stage drugs on target")
    literature_mining_score: Optional[float] = Field(None, description="Score from text-mined literature co-mentions")
    evidence_count: int = Field(default=0, description="Number of evidence datatypes with non-zero scores")
    evidence_breakdown: list[DiseaseLinkEvidence] = Field(
        default_factory=list, description="Per-datatype scores"
    )

    def to_markdown(self) -> str:
        strength = "Strong" if self.overall_score >= 0.5 else "Moderate" if self.overall_score >= 0.2 else "Weak"
        lines = [
            f"## Open Targets: {self.gene_symbol} × {self.disease_name}",
            f"**Overall score:** {self.overall_score:.3f}/1.0 — {strength} evidence ({self.evidence_count} datatypes)",
            f"EFO: `{self.disease_efo_id}` | Ensembl: `{self.ensembl_id}`",
        ]
        if self.evidence_breakdown:
            lines += ["", "| Evidence Type | Score |", "|---|---|"]
            for e in sorted(self.evidence_breakdown, key=lambda x: x.score, reverse=True):
                bar = "█" * int(e.score * 10)
                lines.append(f"| {e.evidence_type} | {e.score:.3f} {bar} |")
        else:
            scores = [
                ("genetic_association", self.genetic_association_score),
                ("somatic_mutation", self.somatic_mutation_score),
                ("known_drug", self.known_drug_score),
                ("literature", self.literature_mining_score),
            ]
            active = [(k, v) for k, v in scores if v is not None]
            if active:
                lines += ["", "| Evidence Type | Score |", "|---|---|"]
                for k, v in sorted(active, key=lambda x: x[1], reverse=True):
                    lines.append(f"| {k} | {v:.3f} |")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# DepMap
# ---------------------------------------------------------------------------


class CellLineEssentiality(BaseModel):
    """CRISPR essentiality score for a single cancer cell line."""

    cell_line: str = Field(description="Cell line name, e.g. 'A375'")
    lineage: str = Field(description="Cancer lineage, e.g. 'Skin'")
    ceres_score: float = Field(description="CERES/Chronos score; < -0.5 indicates dependency")
    is_dependent: bool = Field(description="True if score < -0.5 (standard dependency threshold)")


class CancerDependency(BaseModel):
    """DepMap CRISPR essentiality data for a gene across cancer cell lines."""

    gene_symbol: str
    mean_ceres_score: float = Field(description="Mean score across all profiled cell lines")
    fraction_dependent_lines: float = Field(
        ge=0.0, le=1.0,
        description="Fraction of cell lines with score < -0.5 (dependency threshold)"
    )
    pan_essential: bool = Field(
        description="True if dependent in >90% of lines (core essential — narrow therapeutic window)"
    )
    top_dependent_lineages: list[str] = Field(
        description="Cancer lineages most dependent on this gene (sorted by mean score)"
    )
    cell_lines: list[CellLineEssentiality] = Field(
        default_factory=list, description="Top 10 most dependent cell lines"
    )
    data_source: str = Field(description="DepMap release version or data source description")

    def to_markdown(self) -> str:
        pct = int(self.fraction_dependent_lines * 100)
        if self.pan_essential:
            essentiality_note = f"**PAN-ESSENTIAL** — dependent in {pct}% of lines (narrow therapeutic window)"
        elif pct >= 50:
            essentiality_note = f"**Selective dependency** — {pct}% of cancer lines dependent"
        elif pct >= 20:
            essentiality_note = f"**Partial dependency** — {pct}% of cancer lines dependent"
        else:
            essentiality_note = f"**Weak/no dependency** — {pct}% of cancer lines dependent"

        lines = [
            f"## Cancer Dependency (DepMap): {self.gene_symbol}",
            essentiality_note,
            f"Mean CERES/Chronos score: **{self.mean_ceres_score:.3f}** (threshold: < −0.5 = dependent)",
            f"_Source: {self.data_source}_",
        ]
        if self.top_dependent_lineages:
            lines += ["", f"**Top dependent lineages:** {', '.join(self.top_dependent_lineages[:5])}"]
        if self.cell_lines:
            lines += ["", "| Cell Line | Lineage | Score | Dependent? |", "|---|---|---|---|"]
            for cl in self.cell_lines[:8]:
                dep = "Yes" if cl.is_dependent else "No"
                lines.append(f"| {cl.cell_line} | {cl.lineage} | {cl.ceres_score:.3f} | {dep} |")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# GWAS Catalog
# ---------------------------------------------------------------------------


class GwasHit(BaseModel):
    """A single GWAS association from the GWAS Catalog."""

    study_accession: str = Field(description="GWAS Catalog study accession, e.g. 'GCST000001'")
    trait: str = Field(description="Reported trait or phenotype")
    mapped_gene: str = Field(description="Nearest mapped gene symbol")
    risk_allele: str = Field(description="Risk allele string, e.g. 'rs1234567-A'")
    p_value: float = Field(description="Association p-value (genome-wide significant < 5e-8)")
    beta_or_or: Optional[float] = Field(None, description="Beta coefficient or odds ratio")
    sample_size: Optional[int] = Field(None, description="Total discovery sample size")
    population: Optional[str] = Field(None, description="Ancestry/population of study cohort")
    pubmed_id: Optional[str] = Field(None, description="PubMed ID of the primary publication")


class GwasEvidence(BaseModel):
    """GWAS Catalog associations linking a gene to a trait."""

    gene_symbol: str
    trait_query: str = Field(description="The trait string used for filtering")
    total_associations: int = Field(description="Number of GWAS hits passing the trait filter")
    associations: list[GwasHit] = Field(description="Top associations sorted by p-value")
    strongest_p_value: Optional[float] = Field(None, description="Most significant p-value found")

    def to_markdown(self) -> str:
        if self.total_associations == 0:
            return (
                f"## GWAS Evidence: {self.gene_symbol} × '{self.trait_query}'\n"
                f"No genome-wide significant associations found."
            )
        p_str = f"{self.strongest_p_value:.2e}" if self.strongest_p_value else "N/A"
        lines = [
            f"## GWAS Evidence: {self.gene_symbol} × '{self.trait_query}'",
            f"**{self.total_associations} associations** | Strongest p-value: **{p_str}**",
            "",
            "| Trait | Risk Allele | p-value | Effect Size | N | Study |",
            "|---|---|---|---|---|---|",
        ]
        for hit in self.associations[:8]:
            effect = f"{hit.beta_or_or:.3f}" if hit.beta_or_or is not None else "—"
            n = f"{hit.sample_size:,}" if hit.sample_size else "—"
            lines.append(
                f"| {hit.trait[:40]} | {hit.risk_allele} | {hit.p_value:.2e} | {effect} | {n} | {hit.study_accession} |"
            )
        if len(self.associations) > 8:
            lines.append(f"\n_...and {self.total_associations - 8} more associations_")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# PubChem
# ---------------------------------------------------------------------------


class CompoundActivity(BaseModel):
    """A small molecule with measured bioactivity against the target."""

    cid: int = Field(description="PubChem Compound ID")
    name: str = Field(description="Preferred IUPAC name or common name")
    molecular_formula: Optional[str] = None
    molecular_weight: Optional[float] = Field(None, description="Molecular weight in g/mol")
    activity_outcome: str = Field(description="'Active' | 'Inactive' | 'Inconclusive'")
    activity_value: Optional[float] = Field(None, description="Potency value in nM (IC50/EC50/Ki)")
    activity_type: Optional[str] = Field(None, description="Type of activity measurement, e.g. 'IC50'")
    assay_id: Optional[int] = Field(None, description="PubChem AID of the source assay")


class Compounds(BaseModel):
    """PubChem bioactivity data for small molecules acting on a gene target."""

    gene_symbol: str
    total_active_compounds: int = Field(
        description="Total active compounds found; >50 indicates a well-explored chemical space"
    )
    compounds: list[CompoundActivity] = Field(
        description="Top 20 active compounds sorted by potency (lowest IC50/EC50 first)"
    )

    def to_markdown(self) -> str:
        tractability = "Well-explored" if self.total_active_compounds > 50 else "Emerging"
        lines = [
            f"## PubChem Compounds: {self.gene_symbol}",
            f"**{self.total_active_compounds} active compounds** — {tractability} chemical space",
        ]
        potent = [c for c in self.compounds if c.activity_value is not None]
        if potent:
            lines += ["", "| Compound | MW | Activity | Value (nM) | CID |", "|---|---|---|---|---|"]
            for c in potent[:8]:
                mw = f"{c.molecular_weight:.1f}" if c.molecular_weight else "—"
                val = f"{c.activity_value:.1f}" if c.activity_value else "—"
                atype = c.activity_type or "—"
                lines.append(f"| {c.name[:35]} | {mw} | {atype} | {val} | {c.cid} |")
        elif self.compounds:
            lines += ["", "| Compound | Formula | CID |", "|---|---|---|"]
            for c in self.compounds[:8]:
                formula = c.molecular_formula or "—"
                lines.append(f"| {c.name[:35]} | {formula} | {c.cid} |")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# ChEMBL
# ---------------------------------------------------------------------------


class ChEMBLActivity(BaseModel):
    """A quantitative bioactivity measurement from ChEMBL."""

    molecule_chembl_id: str = Field(description="ChEMBL molecule ID, e.g. 'CHEMBL1'")
    molecule_name: Optional[str] = Field(None, description="Common/trade name if available")
    standard_type: str = Field(description="Measurement type: IC50 | Ki | Kd | EC50")
    pchembl_value: float = Field(
        description="-log10(IC50 in M); 7 = 100 nM lead quality, 9 = 1 nM clinical-grade"
    )
    assay_description: Optional[str] = Field(None, description="Brief assay description from ChEMBL")


class ChEMBLCompounds(BaseModel):
    """ChEMBL bioactivity data for quantitative compound potency against a target."""

    gene_symbol: str
    target_chembl_id: Optional[str] = Field(None, description="ChEMBL target ID, e.g. 'CHEMBL4523582'")
    total_active_compounds: int = Field(
        description="Total compounds with a pChEMBL value for this target"
    )
    best_pchembl: Optional[float] = Field(
        None,
        description="Highest pChEMBL value found (most potent compound); 9 = 1 nM, 7 = 100 nM"
    )
    compounds: list[ChEMBLActivity] = Field(
        description="Top 20 most potent compounds sorted by pChEMBL (descending)"
    )

    def to_markdown(self) -> str:
        if self.best_pchembl is None:
            potency_label = "No quantitative data"
        elif self.best_pchembl >= 9:
            potency_label = f"Clinical-grade (best IC50 ≤ 1 nM, pChEMBL={self.best_pchembl:.1f})"
        elif self.best_pchembl >= 7:
            potency_label = f"Lead quality (best IC50 ≤ 100 nM, pChEMBL={self.best_pchembl:.1f})"
        elif self.best_pchembl >= 5:
            potency_label = f"Hit quality (best IC50 ≤ 10 µM, pChEMBL={self.best_pchembl:.1f})"
        else:
            potency_label = f"Weak (pChEMBL={self.best_pchembl:.1f})"

        lines = [
            f"## ChEMBL Compounds: {self.gene_symbol}",
            f"**{self.total_active_compounds} compounds with potency data** — {potency_label}",
        ]
        if self.target_chembl_id:
            lines.append(f"Target: `{self.target_chembl_id}`")
        if self.compounds:
            lines += ["", "| Compound | Type | pChEMBL | IC50 equiv |", "|---|---|---|---|"]
            for c in self.compounds[:10]:
                name = (c.molecule_name or c.molecule_chembl_id)[:35]
                ic50_nm = 10 ** (9 - c.pchembl_value)
                ic50_str = f"{ic50_nm:.1f} nM" if ic50_nm < 1000 else f"{ic50_nm/1000:.2f} µM"
                lines.append(f"| {name} | {c.standard_type} | {c.pchembl_value:.1f} | {ic50_str} |")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Target prioritization report
# ---------------------------------------------------------------------------


class TargetPrioritizationReport(BaseModel):
    """Structured drug discovery target assessment synthesizing all database evidence."""

    gene_symbol: str
    indication: str
    resolution: GeneResolution
    protein_info: Optional[ProteinInfo] = None
    disease_association: Optional[TargetDiseaseAssociation] = None
    cancer_dependency: Optional[CancerDependency] = None
    gwas_evidence: Optional[GwasEvidence] = None
    compounds: Optional[Compounds] = None
    chembl_compounds: Optional["ChEMBLCompounds"] = None
    priority_score: float = Field(
        ge=0.0, le=10.0,
        description="Composite evidence score 0–10 computed from weighted database evidence"
    )
    priority_tier: str = Field(description="'High' (>7) | 'Medium' (4–7) | 'Low' (<4)")
    evidence_summary: str = Field(
        description="1–3 sentence plain-language summary of key evidence for this target"
    )
    data_gaps: list[str] = Field(
        default_factory=list,
        description="Data sources that returned no data or errors (reduces confidence)"
    )
    errors: dict[str, str] = Field(
        default_factory=dict,
        description="Map of module name -> error message for any failed lookups"
    )

    def to_markdown(self) -> str:
        tier_emoji = {"High": "HIGH PRIORITY", "Medium": "MEDIUM PRIORITY", "Low": "LOW PRIORITY"}
        tier_label = tier_emoji.get(self.priority_tier, self.priority_tier)

        lines = [
            f"# Target Assessment: {self.gene_symbol} | {self.indication}",
            "",
            f"**Priority Score: {self.priority_score:.1f}/10 — {tier_label}**",
            "",
            "## Evidence Summary",
            self.evidence_summary,
            "",
            "## Scoring Breakdown",
            "| Source | Contribution | Max |",
            "|---|---|---|",
        ]

        # Build score breakdown
        da = self.disease_association
        cd = self.cancer_dependency
        gw = self.gwas_evidence
        cp = self.compounds
        pi = self.protein_info

        ot_score = round(da.overall_score * 3.0, 2) if da else 0.0
        lines.append(f"| Open Targets association | {ot_score} | 3.0 |")

        if cd:
            if cd.pan_essential:
                dep_score = 0.5
                dep_note = "pan-essential cap"
            else:
                is_real = "DepMap" in cd.data_source
                confidence = 1.0 if is_real else 0.7
                dep_score = round(cd.fraction_dependent_lines * 2.0 * confidence, 2)
                dep_note = f"{int(cd.fraction_dependent_lines * 100)}% dependent" + ("" if is_real else " (proxy, 0.7×)")
            lines.append(f"| Cancer dependency | {dep_score} ({dep_note}) | 2.0 |")
        else:
            lines.append("| Cancer dependency | 0.0 (no data) | 2.0 |")

        gw_score = round(min(gw.total_associations, 10) / 10 * 2.0, 2) if gw else 0.0
        lines.append(f"| GWAS evidence | {gw_score} | 2.0 |")

        chembl = self.chembl_compounds
        if chembl and chembl.best_pchembl is not None:
            bp = chembl.best_pchembl
            if bp >= 9:
                cp_score = 1.5
            elif bp >= 7:
                cp_score = 1.0
            elif bp >= 5:
                cp_score = 0.5
            else:
                cp_score = 0.25
            cp_label = f"ChEMBL pChEMBL={bp:.1f}"
        elif cp:
            cp_score = round(min(cp.total_active_compounds, 100) / 100 * 1.5, 2)
            cp_label = f"PubChem count={cp.total_active_compounds}"
        else:
            cp_score = 0.0
            cp_label = "no data"
        lines.append(f"| Chemical matter | {cp_score} ({cp_label}) | 1.5 |")

        if pi:
            pi_score = (0.5 if pi.reviewed else 0.0) + min(len(pi.known_variants), 2) / 2 * 1.0
            lines.append(f"| Protein annotation | {pi_score:.2f} | 1.5 |")
        else:
            lines.append("| Protein annotation | 0.0 (no data) | 1.5 |")

        lines += ["", "## Data Sources"]
        sources = [
            ("UniProt", self.protein_info),
            ("Open Targets", self.disease_association),
            ("DepMap", self.cancer_dependency),
            ("GWAS Catalog", self.gwas_evidence),
            ("ChEMBL", self.chembl_compounds),
            ("PubChem", self.compounds),
        ]
        for name, val in sources:
            status = "✓" if val is not None else "✗ (no data)"
            lines.append(f"- **{name}:** {status}")

        if self.data_gaps:
            lines += ["", f"**Data gaps:** {', '.join(self.data_gaps)}"]
        if self.errors:
            lines += ["", "**Errors:**"]
            for module, msg in self.errors.items():
                lines.append(f"- `{module}`: {msg[:120]}")

        # Resolution info
        r = self.resolution
        lines += [
            "",
            "---",
            f"_Resolved: {r.hgnc_symbol}"
            + (f" | NCBI Gene: {r.ncbi_gene_id}" if r.ncbi_gene_id else "")
            + (f" | UniProt: {r.uniprot_accession}" if r.uniprot_accession else "")
            + "_",
        ]

        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Multi-target comparison
# ---------------------------------------------------------------------------


class TargetComparisonRow(BaseModel):
    """Single row in a multi-target comparison."""

    gene_symbol: str
    priority_score: float
    priority_tier: str
    ot_score: Optional[float] = None
    depmap_pct: Optional[int] = None
    depmap_real_data: bool = False
    compound_count: Optional[int] = None
    gwas_count: Optional[int] = None
    data_gaps: list[str] = Field(default_factory=list)
    evidence_summary: str = ""


class ComparisonReport(BaseModel):
    """Side-by-side comparison of multiple drug targets for a given indication."""

    indication: str
    rows: list[TargetComparisonRow]

    def to_markdown(self) -> str:
        ranked = sorted(self.rows, key=lambda r: r.priority_score, reverse=True)
        lines = [
            f"# Target Comparison: {self.indication}",
            "",
            "| Rank | Gene | Score | Tier | OT Score | DepMap % | Compounds | GWAS Hits | Data Gaps |",
            "|---|---|---|---|---|---|---|---|---|",
        ]
        for i, row in enumerate(ranked, 1):
            ot = f"{row.ot_score:.2f}" if row.ot_score is not None else "—"
            dep = f"{row.depmap_pct}%" if row.depmap_pct is not None else "—"
            if row.depmap_pct is not None and not row.depmap_real_data:
                dep += "*"
            cpds = str(row.compound_count) if row.compound_count is not None else "—"
            gwas = str(row.gwas_count) if row.gwas_count is not None else "—"
            gaps = ", ".join(row.data_gaps) if row.data_gaps else "none"
            lines.append(
                f"| {i} | **{row.gene_symbol}** | {row.priority_score:.1f}/10 | {row.priority_tier} "
                f"| {ot} | {dep} | {cpds} | {gwas} | {gaps} |"
            )

        lines += ["", "_* DepMap % estimated from Open Targets somatic mutation data (no direct CRISPR data)_"]
        lines += ["", "## Evidence Summaries"]
        for row in ranked:
            lines += [f"### {row.gene_symbol}", row.evidence_summary, ""]

        return "\n".join(lines)
