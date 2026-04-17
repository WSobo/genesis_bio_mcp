"""Pydantic output models for all MCP tool responses."""

from __future__ import annotations

from pydantic import BaseModel, Field, model_validator

from genesis_bio_mcp.tools.biochem import BiochemFeatures, LiabilityHit

# ---------------------------------------------------------------------------
# Gene resolution
# ---------------------------------------------------------------------------


class GeneResolution(BaseModel):
    """Canonical gene identifiers resolved from a synonym or alias."""

    hgnc_symbol: str = Field(description="Approved HGNC gene symbol, e.g. 'BRAF'")
    hgnc_id: str | None = Field(None, description="HGNC ID, e.g. 'HGNC:1097'")
    ncbi_gene_id: str | None = Field(None, description="NCBI Entrez Gene ID")
    uniprot_accession: str | None = Field(None, description="Primary UniProt accession")
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

    position: str | None = Field(None, description="Sequence position (e.g. '600')")
    original: str | None = Field(None, description="Original amino acid (single letter)")
    variant: str | None = Field(None, description="Variant amino acid (single letter)")
    disease: str | None = Field(None, description="Associated disease name")
    clinical_significance: str | None = Field(
        None, description="Pathogenicity annotation if available"
    )


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
    disulfide_bond_positions: list[int] = Field(
        default_factory=list,
        description="1-indexed Cys positions involved in annotated disulfide bonds (UniProt DISULFID feature).",
    )
    reviewed: bool = Field(description="True if Swiss-Prot (manually reviewed), False if TrEMBL")

    def to_markdown(self) -> str:
        reviewed_label = (
            "Swiss-Prot (manually reviewed)" if self.reviewed else "TrEMBL (unreviewed)"
        )
        lines = [
            f"## Protein: {self.protein_name} ({self.gene_symbol})",
            f"**UniProt:** {self.uniprot_accession} — {reviewed_label}",
            "",
            "### Function",
            self.function_summary,
        ]
        if self.subcellular_locations:
            lines += [
                "",
                f"**Subcellular location:** {', '.join(self.subcellular_locations)}",
            ]
        if self.pathways:
            lines += [
                "",
                f"**Pathways ({len(self.pathways)}):** {', '.join(self.pathways[:5])}",
            ]
            if len(self.pathways) > 5:
                lines[-1] += f" (+{len(self.pathways) - 5} more)"
        if self.disease_associations:
            lines += [
                "",
                f"**Disease associations:** {', '.join(self.disease_associations[:5])}",
            ]
        if self.pdb_structures:
            lines += [
                "",
                f"**PDB structures ({len(self.pdb_structures)}):** {', '.join(self.pdb_structures[:8])}",
            ]
        if self.known_variants:
            lines += ["", "### Disease-linked Variants"]
            for v in self.known_variants[:5]:
                pos = (
                    f"p.{v.original}{v.position}{v.variant}"
                    if (v.original and v.position and v.variant)
                    else v.position or "?"
                )
                disease = v.disease or "unknown disease"
                sig = f" ({v.clinical_significance})" if v.clinical_significance else ""
                lines.append(f"- **{pos}** — {disease}{sig}")
            if len(self.known_variants) > 5:
                lines.append(f"- _...and {len(self.known_variants) - 5} more_")
        if self.disulfide_bond_positions:
            lines += [
                "",
                f"**Annotated disulfide bonds:** {len(self.disulfide_bond_positions) // 2} bond(s) "
                f"across Cys positions {', '.join(str(p) for p in self.disulfide_bond_positions[:12])}"
                + (
                    f" (+{len(self.disulfide_bond_positions) - 12} more)"
                    if len(self.disulfide_bond_positions) > 12
                    else ""
                ),
            ]
        return "\n".join(lines)


class ProteinSequence(BaseModel):
    """Protein sequence with computed biochem features and liability-motif hits."""

    uniprot_accession: str = Field(description="Primary UniProt accession")
    gene_symbol: str = Field(description="HGNC gene symbol")
    organism: str = Field(description="Source organism parsed from FASTA header")
    description: str = Field(description="Free-text description from FASTA header")
    sequence: str = Field(description="Amino-acid sequence (one-letter codes)")
    region_start: int | None = Field(
        None, description="1-indexed slice start, if a region was requested"
    )
    region_end: int | None = Field(
        None, description="1-indexed slice end, if a region was requested"
    )
    features: BiochemFeatures = Field(description="Computed biochemical features")
    liabilities: list[LiabilityHit] = Field(
        default_factory=list, description="Liability-motif scanner hits"
    )
    disulfide_bond_positions: list[int] = Field(
        default_factory=list,
        description="1-indexed Cys positions from UniProt DISULFID annotation",
    )

    def to_markdown(self) -> str:
        region = (
            f" (residues {self.region_start}–{self.region_end})"
            if self.region_start and self.region_end
            else ""
        )
        feats = self.features
        lines = [
            f"## Protein Sequence: {self.gene_symbol}{region}",
            f"**UniProt:** {self.uniprot_accession} | **Organism:** {self.organism or 'n/a'} | "
            f"**Length:** {feats.length} aa",
            "",
            "### Biochemistry",
            "",
            "| Feature | Value |",
            "|---|---|",
            f"| Molecular weight | {feats.molecular_weight_Da:,.1f} Da |",
            f"| Theoretical pI | {feats.theoretical_pI:.2f} |",
            f"| GRAVY (Kyte-Doolittle) | {feats.gravy:+.3f} |",
            f"| Net charge at pH 7.4 | {feats.net_charge_pH74:+.2f} |",
            f"| Aromatic fraction | {feats.aromatic_fraction:.1%} |",
            f"| Cysteines | {feats.cysteine_count} |",
            f"| ε₂₈₀ (reduced) | {feats.extinction_coefficient_280nm_reduced:,} M⁻¹cm⁻¹ |",
            f"| ε₂₈₀ (oxidized, all Cys → S-S) | {feats.extinction_coefficient_280nm:,} M⁻¹cm⁻¹ |",
        ]
        if self.liabilities:
            grouped: dict[str, list[LiabilityHit]] = {}
            for h in self.liabilities:
                grouped.setdefault(h.motif_type, []).append(h)
            lines += ["", "### Liability Motifs"]
            # Order matters: real liabilities first, then observational data.
            order = [
                "deamidation",
                "isomerization",
                "n_glycosylation",
                "free_cysteine",
                "oxidation_methionine",
                "oxidation_tryptophan",
                "cysteine_position",
            ]
            for motif_type in order:
                hits = grouped.get(motif_type)
                if not hits:
                    continue
                positions = ", ".join(str(h.position) for h in hits[:20])
                more = f" (+{len(hits) - 20} more)" if len(hits) > 20 else ""
                label = motif_type.replace("_", " ")
                lines.append(f"- **{label}** ({len(hits)}): positions {positions}{more}")
        if self.disulfide_bond_positions:
            lines += [
                "",
                f"**Annotated disulfide bonds:** {len(self.disulfide_bond_positions) // 2} bond(s) "
                f"across Cys {', '.join(str(p) for p in self.disulfide_bond_positions[:12])}"
                + (
                    f" (+{len(self.disulfide_bond_positions) - 12} more)"
                    if len(self.disulfide_bond_positions) > 12
                    else ""
                ),
            ]
        lines += ["", "### Sequence"]
        if len(self.sequence) <= 80:
            lines.append(f"`{self.sequence}`")
        else:
            lines.append(
                f"`{self.sequence[:60]}...{self.sequence[-20:]}` "
                f"(showing 60 + 20 of {len(self.sequence)} residues)"
            )
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
    overall_score: float = Field(
        ge=0.0, le=1.0, description="Aggregate association score 0–1; >0.5 is strong"
    )
    genetic_association_score: float | None = Field(
        None, description="Score from genetic evidence (GWAS, rare variants)"
    )
    somatic_mutation_score: float | None = Field(
        None, description="Score from somatic mutations in cancer"
    )
    known_drug_score: float | None = Field(
        None, description="Score from approved/clinical-stage drugs on target"
    )
    literature_mining_score: float | None = Field(
        None, description="Score from text-mined literature co-mentions"
    )
    evidence_count: int = Field(
        default=0, description="Number of evidence datatypes with non-zero scores"
    )
    evidence_breakdown: list[DiseaseLinkEvidence] = Field(
        default_factory=list, description="Per-datatype scores"
    )

    def to_markdown(self) -> str:
        strength = (
            "Strong"
            if self.overall_score >= 0.5
            else "Moderate"
            if self.overall_score >= 0.2
            else "Weak"
        )
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
        ge=0.0,
        le=1.0,
        description="Fraction of cell lines with score < -0.5 (dependency threshold)",
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
            essentiality_note = (
                f"**PAN-ESSENTIAL** — dependent in {pct}% of lines (narrow therapeutic window)"
            )
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
            lines += [
                "",
                f"**Top dependent lineages:** {', '.join(self.top_dependent_lineages[:5])}",
            ]
        if self.cell_lines:
            lines += [
                "",
                "| Cell Line | Lineage | Score | Dependent? |",
                "|---|---|---|---|",
            ]
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
    beta_or_or: float | None = Field(None, description="Beta coefficient or odds ratio")
    sample_size: int | None = Field(None, description="Total discovery sample size")
    population: str | None = Field(None, description="Ancestry/population of study cohort")
    pubmed_id: str | None = Field(None, description="PubMed ID of the primary publication")
    efo_uri: str | None = Field(
        None,
        description="EFO ontology URI, e.g. 'http://www.ebi.ac.uk/efo/EFO_0001073' (gene-ID path only)",
    )


class GwasEvidence(BaseModel):
    """GWAS Catalog associations linking a gene to a trait."""

    gene_symbol: str
    trait_query: str = Field(description="The trait string used for filtering")
    total_associations: int = Field(description="Number of GWAS hits passing the trait filter")
    associations: list[GwasHit] = Field(description="Top associations sorted by p-value")
    strongest_p_value: float | None = Field(None, description="Most significant p-value found")

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
    molecular_formula: str | None = None
    molecular_weight: float | None = Field(None, description="Molecular weight in g/mol")
    activity_outcome: str = Field(description="'Active' | 'Inactive' | 'Inconclusive'")
    activity_value: float | None = Field(None, description="Potency value in nM (IC50/EC50/Ki)")
    activity_type: str | None = Field(None, description="Type of activity measurement, e.g. 'IC50'")
    assay_id: int | None = Field(None, description="PubChem AID of the source assay")


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
            lines += [
                "",
                "| Compound | MW | Activity | Value (nM) | CID |",
                "|---|---|---|---|---|",
            ]
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
    molecule_name: str | None = Field(None, description="Common/trade name if available")
    standard_type: str = Field(description="Measurement type: IC50 | Ki | Kd | EC50")
    pchembl_value: float = Field(
        description="-log10(IC50 in M); 7 = 100 nM lead quality, 9 = 1 nM clinical-grade"
    )
    assay_description: str | None = Field(None, description="Brief assay description from ChEMBL")


class ChEMBLCompounds(BaseModel):
    """ChEMBL bioactivity data for quantitative compound potency against a target."""

    gene_symbol: str
    target_chembl_id: str | None = Field(None, description="ChEMBL target ID, e.g. 'CHEMBL4523582'")
    total_active_compounds: int = Field(
        description="Total compounds with a pChEMBL value for this target"
    )
    best_pchembl: float | None = Field(
        None,
        description="Highest pChEMBL value found (most potent compound); 9 = 1 nM, 7 = 100 nM",
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
            lines += [
                "",
                "| Compound | Type | pChEMBL | IC50 equiv |",
                "|---|---|---|---|",
            ]
            for c in self.compounds[:10]:
                name = (c.molecule_name or c.molecule_chembl_id)[:35]
                ic50_nm = 10 ** (9 - c.pchembl_value)
                ic50_str = f"{ic50_nm:.1f} nM" if ic50_nm < 1000 else f"{ic50_nm / 1000:.2f} µM"
                lines.append(f"| {name} | {c.standard_type} | {c.pchembl_value:.1f} | {ic50_str} |")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Protein structure (AlphaFold + RCSB PDB)
# ---------------------------------------------------------------------------


class PDBStructure(BaseModel):
    """A single experimental structure entry from the PDB."""

    pdb_id: str = Field(description="4-character PDB accession, e.g. '4MNE'")
    resolution_angstrom: float | None = Field(None, description="Resolution in Å; lower is better")
    method: str = Field(
        description="Experimental method: X-RAY DIFFRACTION | ELECTRON MICROSCOPY | NMR | etc."
    )
    has_ligand: bool = Field(description="True if a small-molecule ligand is co-crystallized")
    release_year: int | None = Field(None, description="Year the structure was deposited")


class ProteinStructure(BaseModel):
    """Structural data from AlphaFold (predicted) and RCSB PDB (experimental)."""

    gene_symbol: str
    uniprot_accession: str | None = Field(
        None, description="UniProt accession used for AlphaFold lookup"
    )
    alphafold_plddt: float | None = Field(
        None,
        description="Mean pLDDT confidence score 0–100; ≥90 = high confidence, <70 = disordered",
    )
    alphafold_model_url: str | None = Field(None, description="URL to AlphaFold PDB model file")
    alphafold_version: str | None = Field(None, description="AlphaFold model version, e.g. 'v4'")
    experimental_structures: list[PDBStructure] = Field(
        default_factory=list,
        description="Experimental PDB structures sorted by resolution (best first)",
    )
    total_pdb_structures: int = Field(
        default=0, description="Total experimental structures in PDB for this target"
    )
    has_ligand_bound: bool = Field(
        default=False,
        description="True if any PDB structure has a co-crystallized ligand",
    )
    best_resolution: float | None = Field(
        None, description="Best experimental resolution in Å across all structures"
    )

    def to_markdown(self) -> str:
        lines = [f"## Protein Structure: {self.gene_symbol}"]

        # AlphaFold section
        if self.alphafold_plddt is not None:
            conf = (
                "High"
                if self.alphafold_plddt >= 90
                else "Moderate"
                if self.alphafold_plddt >= 70
                else "Low (disordered regions)"
            )
            lines += [
                "",
                f"### AlphaFold Prediction ({self.alphafold_version or 'latest'})",
                f"**Mean pLDDT:** {self.alphafold_plddt:.1f}/100 — {conf} confidence",
            ]
            if self.alphafold_model_url:
                lines.append(f"Model: {self.alphafold_model_url}")
        else:
            lines += ["", "### AlphaFold Prediction", "No AlphaFold model available."]

        # PDB section
        lines += [
            "",
            f"### Experimental Structures (PDB): {self.total_pdb_structures} total",
        ]
        if self.best_resolution:
            lines.append(f"**Best resolution:** {self.best_resolution:.2f} Å")
        if self.has_ligand_bound:
            lines.append(
                "Inhibitor/ligand-bound structures available — structure-based design feasible"
            )
        else:
            lines.append("No ligand-bound structures found — apo structures only")

        if self.experimental_structures:
            lines += [
                "",
                "| PDB ID | Method | Resolution | Ligand Bound | Year |",
                "|---|---|---|---|---|",
            ]
            for s in self.experimental_structures[:8]:
                res = f"{s.resolution_angstrom:.2f} Å" if s.resolution_angstrom else "—"
                lig = "Yes" if s.has_ligand else "No"
                yr = str(s.release_year) if s.release_year else "—"
                lines.append(f"| {s.pdb_id} | {s.method[:25]} | {res} | {lig} | {yr} |")
            if self.total_pdb_structures > 8:
                lines.append(f"\n_...and {self.total_pdb_structures - 8} more structures_")

        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Protein interaction network (STRING)
# ---------------------------------------------------------------------------


class Interactor(BaseModel):
    """A protein interaction partner from STRING."""

    gene_symbol: str = Field(description="HGNC gene symbol of interaction partner")
    protein_name: str | None = Field(None, description="Full protein name")
    score: float = Field(
        ge=0.0,
        le=1.0,
        description="STRING confidence score 0–1; ≥0.9 = high confidence",
    )
    evidence_types: list[str] = Field(
        default_factory=list,
        description="Evidence channels: experiments | database | coexpression | textmining | homology",
    )


class ProteinInteractome(BaseModel):
    """STRING protein interaction network for a gene target."""

    gene_symbol: str
    total_partners: int = Field(description="Total interactors at required_score threshold")
    top_interactors: list[Interactor] = Field(
        description="Top 20 interactors sorted by confidence score (descending)"
    )

    def to_markdown(self) -> str:
        lines = [
            f"## Protein Interactome: {self.gene_symbol}",
            f"**{self.total_partners} interaction partners** (STRING confidence ≥ 0.7)",
        ]
        if self.top_interactors:
            lines += ["", "| Partner | Score | Evidence |", "|---|---|---|"]
            for it in self.top_interactors[:15]:
                ev = ", ".join(it.evidence_types[:4]) if it.evidence_types else "—"
                name = it.protein_name or it.gene_symbol
                lines.append(f"| **{it.gene_symbol}** ({name[:30]}) | {it.score:.3f} | {ev} |")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# BioGRID protein–protein interactions
# ---------------------------------------------------------------------------


class BioGRIDInteraction(BaseModel):
    """A single protein–protein interaction from BioGRID."""

    interactor_a: str = Field(description="Official gene symbol of interactor A")
    interactor_b: str = Field(description="Official gene symbol of interactor B")
    experimental_system: str | None = Field(
        None, description="Experimental method, e.g. 'Two-hybrid', 'Co-IP'"
    )
    experimental_system_type: str | None = Field(
        None, description="Broad type: 'physical' | 'genetic'"
    )
    pubmed_id: str | None = Field(None, description="PubMed ID of the source publication")
    throughput: str | None = Field(
        None, description="Low Throughput | High Throughput | Low/High Throughput"
    )


class BioGRIDInteractome(BaseModel):
    """BioGRID interaction network for a query gene."""

    gene_symbol: str = Field(description="Query gene symbol")
    total_interactions: int = Field(description="Total interactions returned by BioGRID")
    unique_partners: int = Field(description="Number of unique interaction partners")
    interactions: list[BioGRIDInteraction] = Field(
        default_factory=list, description="Individual interaction records (up to 50)"
    )

    def to_markdown(self) -> str:
        lines = [
            f"## BioGRID Interactome: {self.gene_symbol}",
            f"**{self.total_interactions} interactions** | **{self.unique_partners} unique partners**",
        ]
        if self.interactions:
            partner_counts: dict[str, int] = {}
            for ix in self.interactions:
                partner = (
                    ix.interactor_b if ix.interactor_a == self.gene_symbol else ix.interactor_a
                )
                partner_counts[partner] = partner_counts.get(partner, 0) + 1

            lines += [
                "",
                "| Partner | Evidence Count | Method (example) | Type |",
                "|---|---|---|---|",
            ]
            for partner, count in sorted(partner_counts.items(), key=lambda x: -x[1])[:20]:
                example = next(
                    (
                        ix
                        for ix in self.interactions
                        if ix.interactor_a == partner or ix.interactor_b == partner
                    ),
                    None,
                )
                method = example.experimental_system or "—" if example else "—"
                etype = example.experimental_system_type or "—" if example else "—"
                lines.append(f"| **{partner}** | {count} | {method} | {etype} |")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Drug history (DGIdb + ClinicalTrials.gov)
# ---------------------------------------------------------------------------


class DrugInteraction(BaseModel):
    """A known drug-gene interaction from DGIdb."""

    drug_name: str = Field(description="Drug name (common or trade name)")
    interaction_type: str | None = Field(
        None, description="Interaction type: inhibitor | activator | binder | etc."
    )
    phase: int | None = Field(
        None, description="Highest clinical phase: 1 | 2 | 3 | 4 (4 = approved)"
    )
    approved: bool = Field(default=False, description="True if FDA/EMA approved (Phase 4)")
    sources: list[str] = Field(default_factory=list, description="DGIdb data source databases")


class ClinicalTrial(BaseModel):
    """A clinical trial from ClinicalTrials.gov targeting this gene."""

    nct_id: str = Field(description="ClinicalTrials.gov identifier, e.g. 'NCT01234567'")
    title: str = Field(description="Brief trial title")
    phase: str | None = Field(
        None, description="Trial phase: Phase 1 | Phase 2 | Phase 3 | Phase 4 | N/A"
    )
    status: str = Field(
        description="Trial status: RECRUITING | ACTIVE_NOT_RECRUITING | COMPLETED | etc."
    )
    indication: str | None = Field(None, description="Primary condition/indication being studied")


class DrugHistory(BaseModel):
    """Drug development history for a gene target from DGIdb and ClinicalTrials.gov."""

    gene_symbol: str
    known_drugs: list[DrugInteraction] = Field(
        default_factory=list, description="Drugs known to interact with this target"
    )
    approved_drug_count: int = Field(
        default=0, description="Number of FDA/EMA approved drugs on this target"
    )
    trial_counts_by_phase: dict[str, int] = Field(
        default_factory=dict, description="Trial counts keyed by phase string"
    )
    recent_trials: list[ClinicalTrial] = Field(
        default_factory=list, description="Up to 10 most recent trials"
    )

    @model_validator(mode="after")
    def compute_approved_count(self) -> DrugHistory:
        if self.approved_drug_count == 0:
            self.approved_drug_count = sum(1 for d in self.known_drugs if d.approved)
        return self

    def to_markdown(self) -> str:
        total_trials = sum(self.trial_counts_by_phase.values())
        class_label = (
            "First-in-class opportunity"
            if self.approved_drug_count == 0
            else f"Best-in-class ({self.approved_drug_count} approved drugs)"
        )

        lines = [
            f"## Drug History: {self.gene_symbol}",
            f"**{self.approved_drug_count} approved drugs** | {total_trials} clinical trials | {class_label}",
        ]

        if self.known_drugs:
            lines += [
                "",
                "### Known Drugs (DGIdb)",
                "| Drug | Type | Phase | Approved | Sources |",
                "|---|---|---|---|---|",
            ]

            # Sort: approved first, then direct interaction types, then phase desc
            def _drug_sort_key(d: DrugInteraction) -> tuple:
                from genesis_bio_mcp.clients.dgidb import _DIRECT_TYPES

                return (
                    not d.approved,
                    d.interaction_type is None or d.interaction_type.lower() not in _DIRECT_TYPES,
                    -(d.phase or 0),
                    d.drug_name.lower(),
                )

            for d in sorted(self.known_drugs, key=_drug_sort_key)[:12]:
                itype = d.interaction_type or "—"
                phase = f"Phase {d.phase}" if d.phase else "—"
                approved = "Yes" if d.approved else "No"
                sources = ", ".join(d.sources[:3]) if d.sources else "—"
                lines.append(f"| {d.drug_name} | {itype} | {phase} | {approved} | {sources} |")

        if self.trial_counts_by_phase:
            lines += ["", "### Clinical Trial Counts by Phase"]
            for phase, count in sorted(self.trial_counts_by_phase.items()):
                lines.append(f"- {phase}: {count} trial{'s' if count != 1 else ''}")

        if self.recent_trials:
            lines += [
                "",
                "### Recent Trials",
                "| NCT ID | Phase | Status | Indication |",
                "|---|---|---|---|",
            ]
            for t in self.recent_trials[:8]:
                phase = t.phase or "—"
                ind = (t.indication or "—")[:40]
                lines.append(
                    f"| [{t.nct_id}](https://clinicaltrials.gov/study/{t.nct_id}) | {phase} | {t.status} | {ind} |"
                )

        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Pathway context (Reactome)
# ---------------------------------------------------------------------------


class Pathway(BaseModel):
    """A biological pathway from Reactome containing the query gene."""

    reactome_id: str = Field(description="Reactome stable ID, e.g. 'R-HSA-5673001'")
    display_name: str = Field(description="Human-readable pathway name")
    p_value: float | None = Field(None, description="Enrichment p-value (lower = more significant)")
    gene_count: int | None = Field(None, description="Number of genes in this pathway")
    category: str | None = Field(
        None, description="Top-level category: Signaling | Metabolism | Immune | etc."
    )


class PathwayContext(BaseModel):
    """Reactome pathway membership and enrichment for a gene target."""

    gene_symbol: str
    pathways: list[Pathway] = Field(
        description="Pathways sorted by significance (most significant first)"
    )
    top_pathway_name: str | None = Field(None, description="Name of the most significant pathway")

    def to_markdown(self) -> str:
        lines = [
            f"## Pathway Context: {self.gene_symbol}",
            f"**{len(self.pathways)} pathways** found in Reactome",
        ]
        if self.top_pathway_name:
            lines.append(f"Top pathway: **{self.top_pathway_name}**")
        if self.pathways:
            lines += [
                "",
                "| Pathway | Category | Genes | p-value |",
                "|---|---|---|---|",
            ]
            for p in self.pathways[:12]:
                p_str = f"{p.p_value:.2e}" if p.p_value is not None else "—"
                genes = str(p.gene_count) if p.gene_count is not None else "—"
                cat = p.category or "—"
                lines.append(f"| {p.display_name[:50]} | {cat} | {genes} | {p_str} |")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Target prioritization report
# ---------------------------------------------------------------------------


class TargetPrioritizationReport(BaseModel):
    """Structured drug discovery target assessment synthesizing all database evidence."""

    gene_symbol: str
    indication: str
    resolution: GeneResolution
    protein_info: ProteinInfo | None = None
    disease_association: TargetDiseaseAssociation | None = None
    cancer_dependency: CancerDependency | None = None
    gwas_evidence: GwasEvidence | None = None
    compounds: Compounds | None = None
    chembl_compounds: ChEMBLCompounds | None = None
    # Extended mode (populated when extended=True in prioritize_target)
    protein_structure: ProteinStructure | None = None
    protein_interactome: ProteinInteractome | None = None
    drug_history: DrugHistory | None = None
    pathway_context: PathwayContext | None = None
    priority_score: float = Field(
        ge=0.0,
        le=10.0,
        description="Composite evidence score 0–10 computed from weighted database evidence",
    )
    priority_tier: str = Field(description="'High' (>7) | 'Medium' (4–7) | 'Low' (<4)")
    evidence_summary: str = Field(
        description="1–3 sentence plain-language summary of key evidence for this target"
    )
    data_gaps: list[str] = Field(
        default_factory=list,
        description="Data sources that returned no data or errors (reduces confidence)",
    )
    errors: dict[str, str] = Field(
        default_factory=dict,
        description="Map of module name -> error message for any failed lookups",
    )
    # Validation / confidence layer
    data_coverage_pct: float = Field(
        default=0.0,
        ge=0.0,
        le=100.0,
        description="Percentage of the 6 core data sources that returned usable data (0–100)",
    )
    proxy_data_flags: dict[str, bool] = Field(
        default_factory=dict,
        description=(
            "Per-axis flag: True if the score contribution used proxy/estimated data rather than "
            "direct measurement. Keys: 'depmap' (OT somatic proxy instead of CRISPR), "
            "'compounds' (PubChem count only, no ChEMBL potency data)."
        ),
    )
    score_confidence_interval: tuple[float, float] | None = Field(
        default=None,
        description=(
            "(lower_bound, upper_bound) for priority_score based on data completeness. "
            "None when all core sources returned data."
        ),
    )
    api_latency_s: dict[str, float] = Field(
        default_factory=dict,
        description=(
            "Per-API wall-clock latency in seconds for the parallel gather phase. "
            "Keys: 'uniprot', 'open_targets', 'depmap', 'gwas', 'pubchem', 'chembl'."
        ),
    )

    def to_markdown(self) -> str:
        tier_emoji = {
            "High": "HIGH PRIORITY",
            "Medium": "MEDIUM PRIORITY",
            "Low": "LOW PRIORITY",
        }
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

        # Local mirrors of scoring constants from target_prioritization.py — keep in sync.
        _LINEAGE_MATCH_FACTOR = 1.2
        _PUBCHEM_MIN_COMPOUNDS = 5

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
                lineage_match = self.indication and any(
                    self.indication.lower() in lin.lower() or lin.lower() in self.indication.lower()
                    for lin in (cd.top_dependent_lineages or [])
                )
                lineage_factor = _LINEAGE_MATCH_FACTOR if lineage_match else 1.0
                dep_score = round(
                    min(cd.fraction_dependent_lines * 2.0 * confidence * lineage_factor, 2.0), 2
                )
                dep_note = f"{int(cd.fraction_dependent_lines * 100)}% dependent"
                if not is_real:
                    dep_note += " (proxy, 0.7×)"
                if lineage_match:
                    dep_note += f" (lineage match, {_LINEAGE_MATCH_FACTOR}×)"
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
        elif cp and cp.total_active_compounds >= _PUBCHEM_MIN_COMPOUNDS:
            cp_score = round(min(cp.total_active_compounds, 100) / 100 * 1.5, 2)
            cp_label = f"PubChem count={cp.total_active_compounds}"
        elif cp:
            cp_score = 0.0
            cp_label = f"PubChem count={cp.total_active_compounds} (below tractability threshold)"
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

        # Confidence assessment (populated when data_coverage_pct > 0)
        if self.data_coverage_pct > 0 or self.proxy_data_flags or self.score_confidence_interval:
            lines += ["", "## Confidence Assessment"]
            lines.append(
                f"**Data coverage:** {self.data_coverage_pct:.0f}% of core sources returned data"
            )
            if self.score_confidence_interval:
                lo, hi = self.score_confidence_interval
                lines.append(
                    f"**Score range:** {lo:.1f}–{hi:.1f}/10 "
                    f"(uncertainty from {100 - self.data_coverage_pct:.0f}% missing sources)"
                )
            if any(self.proxy_data_flags.values()):
                lines.append("**Proxy data used (reduced confidence):**")
                _proxy_labels = {
                    "depmap": "DepMap — OT somatic mutation proxy instead of real CRISPR data",
                    "compounds": "Chemical matter — PubChem count only, no ChEMBL potency data",
                }
                for axis, is_proxy in self.proxy_data_flags.items():
                    if is_proxy:
                        label = _proxy_labels.get(axis, axis)
                        lines.append(f"  - {label}")

        # API latency breakdown (populated when instrumentation is active)
        if self.api_latency_s:
            slowest = max(self.api_latency_s, key=self.api_latency_s.get)
            lines += ["", "## API Latency"]
            lines.append("| API | Latency (s) |")
            lines.append("|---|---|")
            for api, secs in sorted(self.api_latency_s.items(), key=lambda kv: -kv[1]):
                marker = " ← slowest" if api == slowest else ""
                lines.append(f"| {api} | {secs:.2f}{marker} |")

        # Extended mode sections
        if self.protein_structure:
            lines += ["", self.protein_structure.to_markdown()]
        if self.protein_interactome:
            lines += ["", self.protein_interactome.to_markdown()]
        if self.drug_history:
            lines += ["", self.drug_history.to_markdown()]
        if self.pathway_context:
            lines += ["", self.pathway_context.to_markdown()]

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
    ot_score: float | None = None
    depmap_pct: int | None = None
    depmap_real_data: bool = False
    compound_count: int | None = None
    gwas_count: int | None = None
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

        lines += [
            "",
            "_* DepMap % estimated from Open Targets somatic mutation data (no direct CRISPR data)_",
        ]
        lines += ["", "## Evidence Summaries"]
        for row in ranked:
            lines += [f"### {row.gene_symbol}", row.evidence_summary, ""]

        return "\n".join(lines)


# ---------------------------------------------------------------------------
# IEDB — B-cell epitopes
# ---------------------------------------------------------------------------


class EpitopeRecord(BaseModel):
    """A single B-cell epitope record from IEDB."""

    sequence: str = Field(
        description="Epitope sequence (linear) or discontinuous residue description"
    )
    isotype: str | None = Field(None, description="Antibody isotype, e.g. 'IgG1', 'IgG'")
    pmid: str | None = Field(None, description="PubMed ID of the source publication")
    pdb_id: str | None = Field(None, description="PDB ID if a crystal structure is available")
    antigen_name: str | None = Field(None, description="Antigen protein name from IEDB")
    antigen_accession: str | None = Field(None, description="Antigen UniProt/GenPept accession")
    start_position: int | None = Field(None, description="Epitope start residue position")
    end_position: int | None = Field(None, description="Epitope end residue position")


class EpitopeResults(BaseModel):
    """IEDB B-cell epitope search results for an antigen."""

    antigen_query: str = Field(description="Antigen name used to search IEDB")
    total_assays: int = Field(description="Total positive B-cell assay records retrieved")
    unique_epitopes: int = Field(description="Unique epitope sequences/regions identified")
    with_structure: int = Field(description="Number of assay records with a PDB structure")
    epitopes: list[EpitopeRecord] = Field(
        default_factory=list,
        description="Unique epitope records (up to 50)",
    )

    def to_markdown(self) -> str:
        lines = [
            f"## IEDB B-cell Epitopes — {self.antigen_query}",
            f"**{self.total_assays} positive assays** | "
            f"**{self.unique_epitopes} unique epitopes** | "
            f"{self.with_structure} with PDB structure",
        ]
        if not self.epitopes:
            lines.append("\n_No positive B-cell epitope assays found for this antigen in IEDB._")
            return "\n".join(lines)

        lines += [
            "",
            "| Epitope Sequence / Residues | Isotype | Position | PDB | PubMed |",
            "|---|---|---|---|---|",
        ]
        for ep in self.epitopes[:20]:
            seq = ep.sequence[:60]
            isotype = ep.isotype or "—"
            pos = (
                f"{ep.start_position}–{ep.end_position}"
                if ep.start_position and ep.end_position
                else "—"
            )
            pdb = ep.pdb_id or "—"
            pmid = ep.pmid or "—"
            lines.append(f"| `{seq}` | {isotype} | {pos} | {pdb} | {pmid} |")

        # Key insights
        structural = [ep for ep in self.epitopes if ep.pdb_id]
        linear = [ep for ep in self.epitopes if ep.start_position and ep.end_position]
        insights: list[str] = []
        if structural:
            pdb_ids = ", ".join(set(ep.pdb_id for ep in structural[:5]))  # type: ignore[arg-type]
            insights.append(f"{len(structural)} epitopes with structural evidence (PDB: {pdb_ids})")
        if linear:
            insights.append(f"{len(linear)} epitopes with mapped residue positions")
        if self.unique_epitopes > 10:
            insights.append(
                f"Well-characterized antigen — {self.unique_epitopes} distinct epitopes known"
            )

        if insights:
            lines += ["", "**Key insights:** " + " | ".join(insights)]

        return "\n".join(lines)


# ---------------------------------------------------------------------------
# InterPro — domain annotation
# ---------------------------------------------------------------------------


class DomainAnnotation(BaseModel):
    """A single InterPro domain entry for a protein."""

    interpro_accession: str = Field(description="InterPro accession, e.g. 'IPR000719'")
    name: str = Field(description="Domain name, e.g. 'Protein kinase domain'")
    entry_type: str = Field(
        description="Entry type: domain | family | homologous_superfamily | repeat | site"
    )
    positions: list[tuple[int, int]] = Field(
        default_factory=list,
        description="List of (start, end) residue positions for this domain",
    )
    member_databases: dict[str, list[str]] = Field(
        default_factory=dict,
        description="Member database cross-references, e.g. {'pfam': ['PF00069'], 'smart': ['SM00220']}",
    )
    go_terms: list[str] = Field(
        default_factory=list,
        description="Associated GO term annotations, e.g. 'GO:0004672 protein kinase activity'",
    )


class DomainAnnotations(BaseModel):
    """InterPro domain annotations for a protein."""

    gene_symbol: str = Field(description="HGNC gene symbol")
    uniprot_accession: str = Field(description="UniProt accession used to query InterPro")
    total_entries: int = Field(description="Total InterPro entries for this protein")
    domains: list[DomainAnnotation] = Field(
        default_factory=list,
        description="Domain annotations sorted by start position",
    )

    def to_markdown(self) -> str:
        lines = [
            f"## InterPro Domain Annotations — {self.gene_symbol}",
            f"**UniProt:** {self.uniprot_accession} | **{self.total_entries} InterPro entries**",
        ]
        if not self.domains:
            lines.append("\n_No InterPro domain annotations found for this protein._")
            return "\n".join(lines)

        lines += [
            "",
            "| Domain | Type | Positions | Pfam/SMART | GO Terms |",
            "|---|---|---|---|---|",
        ]
        for d in self.domains:
            pos_str = "; ".join(f"{s}–{e}" for s, e in d.positions[:3])
            pfam = (
                ", ".join(
                    acc
                    for db in ("pfam", "smart", "profile")
                    for acc in d.member_databases.get(db, [])[:2]
                )
                or "—"
            )
            go = "; ".join(t.split(" ", 1)[1] for t in d.go_terms[:2]) if d.go_terms else "—"
            name_short = d.name[:45]
            lines.append(
                f"| **{d.interpro_accession}** {name_short} | {d.entry_type} | {pos_str} | {pfam} | {go} |"
            )

        # Engineering note: list conserved domain regions to avoid
        domain_entries = [d for d in self.domains if d.entry_type == "domain"]
        if domain_entries:
            lines.append("")
            lines.append("**Engineering context (domains to work around):**")
            for d in domain_entries[:5]:
                pos_str = "; ".join(f"aa {s}–{e}" for s, e in d.positions[:2])
                lines.append(f"- **{d.name}** ({pos_str})")

        return "\n".join(lines)


# ---------------------------------------------------------------------------
# gnomAD — evolutionary constraint
# ---------------------------------------------------------------------------


class GnomADConstraint(BaseModel):
    """Gene-level evolutionary constraint metrics from gnomAD v4.

    These metrics quantify how much loss-of-function (LoF) and missense mutation
    the gene tolerates in the human population.  Highly constrained genes (high pLI,
    low oe_lof/LOEUF) are poor candidates for broad mutagenesis; focus engineering
    on tolerant regions or use constraint to guide variant selection.
    """

    gene_symbol: str = Field(description="HGNC gene symbol")
    ensembl_id: str | None = Field(None, description="Ensembl gene ID (GRCh38)")
    gene_name: str | None = Field(None, description="Full gene name")
    constraint_available: bool = Field(
        description="True if gnomAD constraint data exists for this gene"
    )

    # LoF constraint
    pLI: float | None = Field(
        None,
        description=(
            "Probability of loss-of-function intolerance. "
            ">0.9 = highly intolerant (avoid LoF mutations); <0.1 = tolerant."
        ),
    )
    lof_z: float | None = Field(
        None,
        description="LoF Z-score. >3 = significantly constrained against LoF.",
    )
    oe_lof: float | None = Field(
        None,
        description="Observed/expected LoF variant ratio. Lower = more constrained.",
    )
    oe_lof_lower: float | None = Field(
        None, description="Lower bound of oe_lof 90% confidence interval."
    )
    oe_lof_upper: float | None = Field(
        None,
        description=(
            "LOEUF — upper bound of oe_lof confidence interval. "
            "The most widely used constraint metric; <0.35 = highly constrained."
        ),
    )
    obs_lof: int | None = Field(None, description="Observed LoF variants in gnomAD")
    exp_lof: float | None = Field(None, description="Expected LoF variants under neutral model")

    # Missense constraint
    mis_z: float | None = Field(
        None,
        description="Missense Z-score. >3 = significantly depleted of missense variation.",
    )
    oe_mis: float | None = Field(
        None,
        description="Observed/expected missense ratio. <0.6 = significantly constrained.",
    )
    obs_mis: int | None = Field(None, description="Observed missense variants in gnomAD")
    exp_mis: float | None = Field(
        None, description="Expected missense variants under neutral model"
    )

    def to_markdown(self) -> str:
        lines = [f"## gnomAD Constraint — {self.gene_symbol}"]
        if self.gene_name:
            lines[0] += f" ({self.gene_name})"
        if self.ensembl_id:
            lines.append(f"**Ensembl:** {self.ensembl_id}")

        if not self.constraint_available:
            lines.append(
                "\n_No constraint data available in gnomAD for this gene "
                "(insufficient coverage or too few variants)._"
            )
            return "\n".join(lines)

        lines += ["", "### Loss-of-Function Constraint"]
        lines += [
            "| Metric | Value | Interpretation |",
            "|---|---|---|",
        ]
        if self.pLI is not None:
            interp = (
                "Highly intolerant (avoid LoF)"
                if self.pLI > 0.9
                else ("Intermediate" if self.pLI > 0.5 else "Tolerant")
            )
            lines.append(f"| **pLI** | {self.pLI:.4f} | {interp} |")
        if self.oe_lof_upper is not None:
            loeuf_interp = (
                "Highly constrained"
                if self.oe_lof_upper < 0.35
                else ("Constrained" if self.oe_lof_upper < 0.6 else "Unconstrained")
            )
            lines.append(f"| **LOEUF** (oe_lof upper) | {self.oe_lof_upper:.3f} | {loeuf_interp} |")
        if self.oe_lof is not None:
            lines.append(f"| oe_lof | {self.oe_lof:.3f} | — |")
        if self.lof_z is not None:
            lines.append(
                f"| LoF Z-score | {self.lof_z:.2f} | {'>3 = constrained' if self.lof_z > 3 else 'Normal'} |"
            )
        if self.obs_lof is not None and self.exp_lof is not None:
            lines.append(f"| obs/exp LoF | {self.obs_lof} / {self.exp_lof:.1f} | — |")

        lines += ["", "### Missense Constraint"]
        lines += [
            "| Metric | Value | Interpretation |",
            "|---|---|---|",
        ]
        if self.oe_mis is not None:
            mis_interp = "Constrained" if self.oe_mis < 0.6 else "Tolerant"
            lines.append(f"| **oe_mis** | {self.oe_mis:.3f} | {mis_interp} |")
        if self.mis_z is not None:
            lines.append(
                f"| Missense Z-score | {self.mis_z:.2f} | {'>3 = constrained' if self.mis_z > 3 else 'Normal'} |"
            )
        if self.obs_mis is not None and self.exp_mis is not None:
            lines.append(f"| obs/exp missense | {self.obs_mis} / {self.exp_mis:.1f} | — |")

        # Engineering implication
        lines.append("")
        if self.pLI is not None and self.oe_lof_upper is not None:
            if self.pLI > 0.9 or self.oe_lof_upper < 0.35:
                lines.append(
                    "**Engineering note:** High constraint — this gene does not tolerate LoF "
                    "mutations in the population. Focus on tolerant missense positions and avoid "
                    "disrupting essential domains."
                )
            elif self.oe_mis is not None and self.oe_mis < 0.6:
                lines.append(
                    "**Engineering note:** Moderate missense constraint. "
                    "Prefer conservative substitutions; validate fitness experimentally."
                )
            else:
                lines.append(
                    "**Engineering note:** Relatively unconstrained — tolerates both LoF and "
                    "missense variation. Broader mutagenesis strategies are feasible."
                )

        return "\n".join(lines)


# ---------------------------------------------------------------------------
# SAbDab — antibody / nanobody structures
# ---------------------------------------------------------------------------


class AntibodyStructure(BaseModel):
    """A single antibody or nanobody (VHH) structure from SAbDab."""

    pdb: str = Field(description="PDB accession, e.g. '7KMG'")
    is_nanobody: bool = Field(
        description="True if this is a VHH/nanobody (single-domain, no light chain)"
    )
    antigen_name: str | None = Field(None, description="Antigen name as annotated in SAbDab")
    resolution_ang: float | None = Field(
        None, description="Crystal/cryo-EM resolution in Ångströms (lower = better)"
    )
    method: str | None = Field(
        None, description="Experimental method: X-RAY DIFFRACTION | ELECTRON MICROSCOPY | NMR"
    )
    heavy_species: str | None = Field(
        None, description="Species of the heavy chain (antibody organism)"
    )
    light_species: str | None = Field(
        None, description="Species of the light chain; None for nanobodies"
    )
    heavy_subclass: str | None = Field(None, description="IGHV germline family, e.g. 'IGHV3'")
    light_subclass: str | None = Field(
        None, description="IGKV/IGLV germline family; None for nanobodies"
    )
    is_engineered: bool = Field(
        description="True if the antibody sequence is engineered / recombinant"
    )
    is_scfv: bool = Field(description="True if this is a single-chain Fv fragment")
    affinity_nM: float | None = Field(
        None, description="Binding affinity in nM if experimentally measured"
    )
    compound: str | None = Field(None, description="PDB compound description")
    date_added: str | None = Field(None, description="Date added to SAbDab (MM/DD/YY)")
    pmid: str | None = Field(None, description="PubMed ID of the source publication")
    vh_cdr1: str | None = Field(None, description="VH CDR1 sequence (Chothia)")
    vh_cdr2: str | None = Field(None, description="VH CDR2 sequence (Chothia)")
    vh_cdr3: str | None = Field(None, description="VH CDR3 sequence (Chothia)")
    vl_cdr1: str | None = Field(None, description="VL CDR1 sequence (Chothia); None for nanobodies")
    vl_cdr2: str | None = Field(None, description="VL CDR2 sequence (Chothia); None for nanobodies")
    vl_cdr3: str | None = Field(None, description="VL CDR3 sequence (Chothia); None for nanobodies")


class AntibodyStructures(BaseModel):
    """SAbDab search results: antibody and nanobody structures for an antigen."""

    query: str = Field(description="Antigen query string used to search SAbDab")
    total_structures: int = Field(description="Total structures matching the query")
    nanobody_count: int = Field(description="Number of VHH/nanobody entries in the results")
    fab_count: int = Field(description="Number of conventional Fab/IgG entries in the results")
    structures: list[AntibodyStructure] = Field(
        default_factory=list,
        description="Top structures sorted by resolution (best first)",
    )

    def to_markdown(self) -> str:
        lines = [
            f"## Antibody Structures — {self.query}",
            f"**{self.total_structures} structures** | "
            f"{self.fab_count} Fab/IgG | {self.nanobody_count} nanobody (VHH)",
        ]
        if not self.structures:
            lines.append("\n_No antibody or nanobody structures found in SAbDab for this query._")
            return "\n".join(lines)

        lines += [
            "",
            "| PDB | Type | Resolution | Method | Ab Species | Subclass | Engineered | Affinity (nM) |",
            "|---|---|---|---|---|---|---|---|",
        ]
        for s in self.structures:
            ab_type = "VHH" if s.is_nanobody else ("scFv" if s.is_scfv else "Fab/IgG")
            res = f"{s.resolution_ang:.2f} Å" if s.resolution_ang else "—"
            method = (s.method or "—")[:12]
            species = (s.heavy_species or "—")[:20]
            subclass = s.heavy_subclass or "—"
            eng = "Yes" if s.is_engineered else "No"
            aff = f"{s.affinity_nM:.1f}" if s.affinity_nM else "—"
            lines.append(
                f"| **{s.pdb}** | {ab_type} | {res} | {method} | {species} | {subclass} | {eng} | {aff} |"
            )

        # CDR annotation block for top structure
        top = self.structures[0]
        if top.vh_cdr1 or top.vh_cdr2 or top.vh_cdr3:
            lines += ["", f"**CDR sequences ({top.pdb}, Chothia):**"]
            if top.vh_cdr1 or top.vh_cdr2 or top.vh_cdr3:
                lines.append(
                    f"VH: CDR1=`{top.vh_cdr1 or '—'}` | CDR2=`{top.vh_cdr2 or '—'}` | CDR3=`{top.vh_cdr3 or '—'}`"
                )
            if not top.is_nanobody and (top.vl_cdr1 or top.vl_cdr2 or top.vl_cdr3):
                lines.append(
                    f"VL: CDR1=`{top.vl_cdr1 or '—'}` | CDR2=`{top.vl_cdr2 or '—'}` | CDR3=`{top.vl_cdr3 or '—'}`"
                )

        # Surface key insights
        best = self.structures[0]
        insights: list[str] = []
        if best.resolution_ang:
            insights.append(f"Best resolution: **{best.pdb}** at {best.resolution_ang:.2f} Å")
        if self.nanobody_count > 0:
            insights.append(
                f"{self.nanobody_count} VHH nanobod{'ies' if self.nanobody_count > 1 else 'y'} available"
            )
        affinities = [s.affinity_nM for s in self.structures if s.affinity_nM is not None]
        if affinities:
            best_aff = min(affinities)
            insights.append(f"Best measured affinity: {best_aff:.1f} nM")

        if insights:
            lines += ["", "**Key insights:** " + " | ".join(insights)]

        return "\n".join(lines)


# ---------------------------------------------------------------------------
# MaveDB — deep mutational scanning score sets
# ---------------------------------------------------------------------------


class DMSScoreSet(BaseModel):
    """A single deep mutational scanning (DMS) score set from MaveDB."""

    urn: str = Field(description="MaveDB URN identifier, e.g. 'urn:mavedb:00000001-a-1'")
    title: str = Field(description="Score set title describing the experiment")
    short_description: str | None = Field(None, description="Brief description of the assay")
    num_variants: int = Field(description="Number of sequence variants scored in this set")
    target_gene: str | None = Field(None, description="Target gene symbol from MaveDB metadata")
    uniprot_accession: str | None = Field(
        None, description="UniProt accession of the target protein if mapped"
    )
    published_date: str | None = Field(None, description="Publication date of the score set")
    pmid: str | None = Field(None, description="PubMed ID of the primary publication")
    doi: str | None = Field(None, description="DOI of the primary publication")


class DMSResults(BaseModel):
    """Deep mutational scanning datasets from MaveDB for a gene."""

    gene_symbol: str = Field(description="Query gene symbol")
    total_score_sets: int = Field(description="Total number of DMS score sets found")
    total_variants: int = Field(description="Sum of variants across all score sets")
    score_sets: list[DMSScoreSet] = Field(
        default_factory=list,
        description="Available DMS score sets sorted by variant count (largest first)",
    )

    def to_markdown(self) -> str:
        lines = [
            f"## Deep Mutational Scanning — {self.gene_symbol}",
            f"**{self.total_score_sets} score set(s)** | "
            f"**{self.total_variants:,} total variants** across all experiments",
        ]
        if not self.score_sets:
            lines.append(
                "\n_No DMS score sets found in MaveDB for this gene. "
                "DMS data is sparse — not all genes have been characterized._"
            )
            return "\n".join(lines)

        lines += [
            "",
            "| URN | Title | Variants | UniProt | Published | PMID |",
            "|---|---|---|---|---|---|",
        ]
        for ss in self.score_sets:
            title = ss.title[:55] + "…" if len(ss.title) > 55 else ss.title
            uniprot = ss.uniprot_accession or "—"
            pub = ss.published_date or "—"
            pmid = ss.pmid or "—"
            lines.append(
                f"| `{ss.urn}` | {title} | {ss.num_variants:,} | {uniprot} | {pub} | {pmid} |"
            )

        lines += [
            "",
            "**Note:** Retrieve individual variant scores via the MaveDB web interface "
            "or API using the URN. Higher variant counts indicate more complete mutational coverage.",
        ]
        return "\n".join(lines)
