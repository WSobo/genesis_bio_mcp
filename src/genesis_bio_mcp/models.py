"""Pydantic output models for all MCP tool responses."""

from __future__ import annotations

from typing import Optional
from pydantic import BaseModel, Field


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
