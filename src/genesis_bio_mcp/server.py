"""genesis_bio_mcp MCP server.

Exposes 15 tools for biomedical database queries:
  - resolve_gene                  UniProt + NCBI: gene symbol → canonical IDs
  - get_protein_info              UniProt Swiss-Prot protein annotation
  - get_target_disease_association Open Targets: target–disease association score
  - get_cancer_dependency         DepMap: CRISPR essentiality across cancer lines
  - get_gwas_evidence             GWAS Catalog: genetic associations for a trait
  - get_compounds                 PubChem: active small molecules against a target
  - get_chembl_compounds          ChEMBL: quantitative IC50/Ki/Kd potency data
  - get_protein_structure         AlphaFold + RCSB PDB: structural data
  - get_protein_interactome       STRING: binding partners and selectivity risks
  - get_drug_history              DGIdb + ClinicalTrials.gov: known drugs and trials
  - get_pathway_context           Reactome: pathway membership and enrichment for a gene
  - get_pathway_members           Reactome: enumerate all genes in a named pathway
  - prioritize_target             Orchestration: full target assessment report
  - compare_targets               Compare 2–5 targets side by side for an indication
  - run_biology_workflow          AI agent: dynamic tool selection for multi-step questions

All tools return Markdown strings for direct LLM consumption.
"""

from __future__ import annotations

import asyncio
import json as _json
import logging
import os
from contextlib import asynccontextmanager
from types import SimpleNamespace
from typing import Literal

import httpx
from mcp.server.fastmcp import FastMCP
from mcp.types import ToolAnnotations
from pydantic import BaseModel, ConfigDict, Field

from genesis_bio_mcp import __version__
from genesis_bio_mcp.clients.alphafold import AlphaFoldClient
from genesis_bio_mcp.clients.chembl import ChEMBLClient
from genesis_bio_mcp.clients.clinical_trials import ClinicalTrialsClient
from genesis_bio_mcp.clients.depmap import DepMapClient, load_depmap_cache
from genesis_bio_mcp.clients.dgidb import DGIdbClient
from genesis_bio_mcp.clients.gwas import GwasClient
from genesis_bio_mcp.clients.open_targets import OpenTargetsClient
from genesis_bio_mcp.clients.pubchem import PubChemClient
from genesis_bio_mcp.clients.reactome import ReactomeClient
from genesis_bio_mcp.clients.string_db import StringDbClient
from genesis_bio_mcp.clients.uniprot import UniProtClient
from genesis_bio_mcp.config.efo_resolver import EFOResolver
from genesis_bio_mcp.models import ComparisonReport, DrugHistory, TargetComparisonRow
from genesis_bio_mcp.tools.gene_resolver import resolve_gene as _resolve_gene
from genesis_bio_mcp.tools.target_prioritization import (
    prioritize_target as _prioritize_target,
)
from genesis_bio_mcp.workflow_agent import (
    build_tool_registry,
    format_registry_docs,
    run_agent_loop,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

_HEADERS = {
    "User-Agent": f"genesis-bio-mcp/{__version__} (research; github.com/WSobo/genesis-bio-mcp)",
    "Accept": "application/json",
}


@asynccontextmanager
async def lifespan(server: FastMCP):
    """Manage a shared httpx.AsyncClient and pre-load the DepMap gene cache."""
    if not os.environ.get("ANTHROPIC_API_KEY"):
        logger.warning(
            "ANTHROPIC_API_KEY is not set. run_biology_workflow will fail. "
            "Set it in claude_desktop_config.json under 'env' or export it before "
            "starting the MCP server."
        )

    async with httpx.AsyncClient(headers=_HEADERS, timeout=30.0, follow_redirects=True) as client:
        # Fetch DepMap gene_dep_summary once at startup for instant lookups
        gene_dep_cache = await load_depmap_cache(client)

        server.state = SimpleNamespace()
        server.state.uniprot = UniProtClient(client)
        server.state.open_targets = OpenTargetsClient(client)
        server.state.depmap = DepMapClient(client, gene_dep_cache)
        server.state.gwas = GwasClient(client, efo_resolver=EFOResolver(client))
        server.state.pubchem = PubChemClient(client)
        server.state.chembl = ChEMBLClient(client)
        server.state.alphafold = AlphaFoldClient(client)
        server.state.string_db = StringDbClient(client)
        server.state.dgidb = DGIdbClient(client)
        server.state.clinical_trials = ClinicalTrialsClient(client)
        server.state.reactome = ReactomeClient(client)
        yield


mcp = FastMCP("genesis_bio_mcp", lifespan=lifespan)


# ---------------------------------------------------------------------------
# Internal helper: alias-tolerant symbol resolution
# ---------------------------------------------------------------------------


async def _resolve_symbol(gene_name: str) -> tuple[str, str | None]:
    """Resolve a gene name or alias to (canonical_hgnc_symbol, ncbi_gene_id).

    Called at the top of every individual tool so that common aliases
    (HER2 → ERBB2, p53 → TP53, COX2 → PTGS2) are transparently resolved
    before any database query.  Falls back silently to the uppercased input
    if resolution fails, so tools never hard-fail on lookup errors.

    Returns:
        (hgnc_symbol, ncbi_gene_id) — ncbi_gene_id may be None if NCBI lookup failed.
    """
    try:
        resolution = await _resolve_gene(gene_name, uniprot_client=mcp.state.uniprot)
        symbol = resolution.hgnc_symbol or gene_name.strip().upper()
        return symbol, resolution.ncbi_gene_id
    except Exception:
        return gene_name.strip().upper(), None


def _fmt(result: object, response_format: str, error_msg: str) -> str:
    """Format a Pydantic model as Markdown or JSON, or return an error representation.

    Args:
        result: Pydantic model with .to_markdown() and .model_dump_json(), or None.
        response_format: "markdown" (default) or "json".
        error_msg: Human-readable error used when result is None.
    """
    if result is None:
        if response_format == "json":
            return _json.dumps({"error": error_msg})
        return f"**Error:** {error_msg}"
    if response_format == "json":
        return result.model_dump_json(indent=2)  # type: ignore[attr-defined]
    return result.to_markdown()  # type: ignore[attr-defined]


# ---------------------------------------------------------------------------
# Input models — Pydantic V2 with strict validation and auto whitespace strip
# ---------------------------------------------------------------------------

_RESPONSE_FORMAT_FIELD = Field(
    default="markdown",
    description="Output format: 'markdown' (human-readable) or 'json' (machine-readable model).",
)
_GENE_SYMBOL_FIELD = Field(
    ...,
    description="HGNC gene symbol. Examples: 'BRAF', 'EGFR', 'TP53', 'KRAS', 'PCSK9'.",
    min_length=1,
    max_length=50,
)


class _GeneInput(BaseModel):
    """Shared base for single-gene lookup tools."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    gene_symbol: str = _GENE_SYMBOL_FIELD
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD


class GetProteinInfoInput(_GeneInput):
    """Input for get_protein_info."""


class GetCancerDependencyInput(_GeneInput):
    """Input for get_cancer_dependency."""


class GetCompoundsInput(_GeneInput):
    """Input for get_compounds."""


class GetChEMBLCompoundsInput(_GeneInput):
    """Input for get_chembl_compounds."""


class GetProteinStructureInput(_GeneInput):
    """Input for get_protein_structure."""


class GetProteinInteractomeInput(_GeneInput):
    """Input for get_protein_interactome."""


class GetDrugHistoryInput(_GeneInput):
    """Input for get_drug_history."""


class GetPathwayContextInput(_GeneInput):
    """Input for get_pathway_context."""


class GetPathwayMembersInput(BaseModel):
    """Input for get_pathway_members."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    pathway_name_or_id: str = Field(
        ...,
        description=(
            "Reactome pathway display name (e.g. 'MAPK signaling', 'RAF/MAP kinase cascade') "
            "or stable ID (e.g. 'R-HSA-5673001')."
        ),
        min_length=1,
        max_length=200,
    )
    max_genes: int = Field(
        default=50,
        description="Maximum number of gene symbols to return (default 50).",
        ge=1,
        le=500,
    )


class ResolveGeneInput(BaseModel):
    """Input for resolve_gene — uses gene_name (accepts aliases) rather than a canonical symbol."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    gene_name: str = Field(
        ...,
        description=(
            "Gene symbol, alias, or synonym (case-insensitive). "
            "Examples: 'BRAF', 'HER2', 'p53', 'B-RAF', 'ErbB2'."
        ),
        min_length=1,
        max_length=50,
    )
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD


class GetTargetDiseaseInput(BaseModel):
    """Input for get_target_disease_association."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    gene_symbol: str = _GENE_SYMBOL_FIELD
    disease_name: str = Field(
        ...,
        description=(
            "Free-text disease or indication name. Open Targets maps this to the closest "
            "EFO ontology term. Examples: 'melanoma', 'non-small cell lung cancer', "
            "'type 2 diabetes', 'Alzheimer disease'."
        ),
        min_length=1,
        max_length=200,
    )
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD


class GetGwasEvidenceInput(BaseModel):
    """Input for get_gwas_evidence."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    gene_symbol: str = _GENE_SYMBOL_FIELD
    trait: str = Field(
        ...,
        description=(
            "Trait or disease name for filtering GWAS associations (case-insensitive substring "
            "match). Examples: 'type 2 diabetes', 'body mass index', 'breast cancer'."
        ),
        min_length=1,
        max_length=200,
    )
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD


class PrioritizeTargetInput(BaseModel):
    """Input for prioritize_target."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    gene_symbol: str = _GENE_SYMBOL_FIELD
    indication: str = Field(
        ...,
        description=(
            "Therapeutic indication or disease area. "
            "Examples: 'melanoma', 'non-small cell lung cancer', 'type 2 diabetes'."
        ),
        min_length=1,
        max_length=200,
    )
    extended: bool = Field(
        default=False,
        description=(
            "If True, also retrieves AlphaFold structure, STRING interactome, drug history, "
            "and Reactome pathway context. Adds ~10–20 s to query time."
        ),
    )
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD


class CompareTargetsInput(BaseModel):
    """Input for compare_targets."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    gene_symbols: list[str] = Field(
        ...,
        description="List of 2–5 HGNC gene symbols. Example: ['BRAF', 'EGFR', 'KRAS'].",
        min_length=2,
    )
    indication: str = Field(
        ...,
        description=(
            "Therapeutic indication shared across all targets. "
            "Example: 'melanoma', 'non-small cell lung cancer'."
        ),
        min_length=1,
        max_length=200,
    )
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD


class RunBiologyWorkflowInput(BaseModel):
    """Input for run_biology_workflow."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    question: str = Field(
        ...,
        description=(
            "Free-text biology or drug discovery question requiring multi-step reasoning. "
            "Examples: 'Find underexplored targets in the MAPK pathway with no approved drugs', "
            "'Is KRAS druggable? What is the best chemical matter available?'."
        ),
        min_length=5,
        max_length=2000,
    )


# ---------------------------------------------------------------------------
# Tool definitions — all return Markdown strings for LLM readability
# ---------------------------------------------------------------------------


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def resolve_gene(params: ResolveGeneInput) -> str:
    """Resolve a gene name or alias to canonical identifiers across databases.

    Use this tool FIRST when the input is a gene alias, synonym, or informal name.
    Provides HGNC symbol, NCBI Gene ID, and UniProt accession needed by other tools.

    Args:
        params (ResolveGeneInput): gene_name, response_format.

    Returns:
        Markdown with canonical symbol, NCBI Gene ID, UniProt accession, and synonyms.
    """
    result = await _resolve_gene(params.gene_name, uniprot_client=mcp.state.uniprot)
    return _fmt(result, params.response_format, f"Could not resolve gene '{params.gene_name}'")


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_protein_info(params: GetProteinInfoInput) -> str:
    """Retrieve protein-level annotation for a human gene from UniProt Swiss-Prot.

    Use this tool to understand a protein's biological function, subcellular location,
    known disease-linked variants, and available 3D structures. Best used after
    resolve_gene to ensure a canonical symbol.

    Args:
        params (GetProteinInfoInput): gene_symbol, response_format.

    Returns:
        Markdown with function summary, pathways, disease associations, PDB IDs,
        known variants, and reviewed status.
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    result = await mcp.state.uniprot.get_protein(symbol)
    return _fmt(
        result,
        params.response_format,
        f"No UniProt Swiss-Prot entry found for gene '{symbol}' in Homo sapiens.",
    )


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_target_disease_association(params: GetTargetDiseaseInput) -> str:
    """Query Open Targets for the evidence-based association score between a gene and disease.

    Use this tool to assess genetic, clinical, and literature evidence linking a drug
    target to a specific indication. The overall_score (0–1) is Open Targets' aggregate
    evidence strength; scores >0.5 are considered strong support for a target–disease link.

    Args:
        params (GetTargetDiseaseInput): gene_symbol, disease_name, response_format.

    Returns:
        Markdown with overall_score and per-datatype evidence scores
        (genetic_association, somatic_mutation, known_drug, literature).
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    result = await mcp.state.open_targets.get_association(symbol, params.disease_name)
    return _fmt(
        result,
        params.response_format,
        f"No Open Targets association found for '{symbol}' and '{params.disease_name}'.",
    )


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_cancer_dependency(params: GetCancerDependencyInput) -> str:
    """Retrieve CRISPR essentiality scores for a gene across cancer cell lines from DepMap.

    Uses real DepMap Chronos Combined data when available (loaded at server startup),
    supplemented by Open Targets somatic mutation data for lineage context.
    Pan-essential genes (common_essential=True) are core cellular machinery and may
    have narrow therapeutic windows.

    Args:
        params (GetCancerDependencyInput): gene_symbol, response_format.

    Returns:
        Markdown with fraction of dependent lines, pan-essential flag, top lineages,
        and the data source (real DepMap or OT proxy).
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    result = await mcp.state.depmap.get_essentiality(symbol)
    return _fmt(
        result,
        params.response_format,
        f"DepMap essentiality data unavailable for '{symbol}'.",
    )


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_gwas_evidence(params: GetGwasEvidenceInput) -> str:
    """Retrieve GWAS Catalog associations linking a gene locus to a disease trait.

    Use this tool to find genome-wide significant SNP associations (p < 5e-8) near
    a gene for a phenotypic trait. High association counts and low p-values strengthen
    genetic causality evidence for target selection.

    Args:
        params (GetGwasEvidenceInput): gene_symbol, trait, response_format.

    Returns:
        Markdown with GWAS hit count, strongest p-value, and top associations table.
    """
    symbol, ncbi_gene_id = await _resolve_symbol(params.gene_symbol)
    result = await mcp.state.gwas.get_evidence(symbol, params.trait, ncbi_gene_id=ncbi_gene_id)
    return _fmt(
        result,
        params.response_format,
        f"No GWAS Catalog associations found for gene '{symbol}' and trait '{params.trait}'.",
    )


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_compounds(params: GetCompoundsInput) -> str:
    """Retrieve small molecules with bioactivity against a gene target from PubChem.

    Use this tool to assess target druggability — whether active chemical matter exists.
    Returns active compounds sorted by potency (lowest IC50/EC50 first). A count >50
    indicates a well-explored chemical space and tractable target.

    Args:
        params (GetCompoundsInput): gene_symbol, response_format.

    Returns:
        Markdown with total active compound count and top compounds by potency (IC50/EC50 in nM).
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    result = await mcp.state.pubchem.get_compounds(symbol)
    return _fmt(
        result, params.response_format, f"No PubChem bioactivity data found for gene '{symbol}'."
    )


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_chembl_compounds(params: GetChEMBLCompoundsInput) -> str:
    """Retrieve quantitative potency data (IC50/Ki/Kd) from ChEMBL for a gene target.

    Use this tool when you need potency-based compound data rather than just activity
    counts. ChEMBL's pChEMBL values (−log10 of molar IC50/Ki/Kd) allow direct
    compound ranking: pChEMBL ≥ 9 = clinical-grade (≤1 nM), ≥ 7 = lead quality,
    ≥ 5 = hit quality. Complements get_compounds (PubChem) which reports activity
    count but lacks standardized potency metrics.

    Args:
        params (GetChEMBLCompoundsInput): gene_symbol, response_format.

    Returns:
        Markdown with ChEMBL target ID, best pChEMBL value, compound count, and
        a ranked table of top compounds with assay types and potency values.
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    result = await mcp.state.chembl.get_compounds(symbol)
    return _fmt(
        result, params.response_format, f"No ChEMBL bioactivity data found for gene '{symbol}'."
    )


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_protein_structure(params: GetProteinStructureInput) -> str:
    """Retrieve structural data for a protein from AlphaFold and RCSB PDB.

    Use this tool to assess structural feasibility for drug design. Reports AlphaFold
    prediction confidence (pLDDT), experimental PDB structures, best resolution,
    and whether inhibitor-bound structures exist (critical for structure-based design).

    Args:
        params (GetProteinStructureInput): gene_symbol, response_format.

    Returns:
        Markdown with AlphaFold pLDDT score, PDB structure count and resolution,
        and whether ligand-bound structures are available.
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    protein = await mcp.state.uniprot.get_protein(symbol)
    accession = protein.uniprot_accession if protein else None
    result = await mcp.state.alphafold.get_structure(symbol, uniprot_accession=accession)
    return _fmt(
        result,
        params.response_format,
        f"No structural data found for '{symbol}'. Ensure the gene symbol is correct.",
    )


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_protein_interactome(params: GetProteinInteractomeInput) -> str:
    """Retrieve protein interaction partners from STRING to assess selectivity risks.

    Use this tool to understand which proteins interact with your target — binding
    partners, paralogs, and pathway neighbors that your compound might also engage.
    High-confidence interactors (score ≥ 0.9) are particularly important for selectivity.

    Args:
        params (GetProteinInteractomeInput): gene_symbol, response_format.

    Returns:
        Markdown with top 20 interaction partners sorted by STRING confidence score,
        evidence types (experiments, database, coexpression), and total partner count.
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    result = await mcp.state.string_db.get_interactome(symbol)
    if result is not None and result.total_partners == 0:
        result = None
    return _fmt(result, params.response_format, f"No STRING interaction data found for '{symbol}'.")


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_drug_history(params: GetDrugHistoryInput) -> str:
    """Retrieve the drug development history for a gene target.

    Use this tool to understand the clinical precedent for targeting a gene:
    what drugs already exist (approved or investigational), what indications have
    been pursued, and how many clinical trials have targeted this gene. Essential
    for first-in-class vs. best-in-class strategy decisions.

    Args:
        params (GetDrugHistoryInput): gene_symbol, response_format.

    Returns:
        Markdown with known drugs (type, phase, approval status), trial counts by
        phase, and a table of recent clinical trials from ClinicalTrials.gov.
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    drugs, ct_result = await asyncio.gather(
        mcp.state.dgidb.get_drug_interactions(symbol),
        mcp.state.clinical_trials.get_trials(symbol),
    )
    ct_trials, ct_counts = ct_result
    result = DrugHistory(
        gene_symbol=symbol,
        known_drugs=drugs,
        approved_drug_count=sum(1 for d in drugs if d.approved),
        trial_counts_by_phase=ct_counts,
        recent_trials=ct_trials[:10],
    )
    if params.response_format == "json":
        return result.model_dump_json(indent=2)
    if not drugs and not ct_trials:
        return (
            f"**No drug history found for '{symbol}'.** This may be a first-in-class opportunity."
        )
    return result.to_markdown()


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_pathway_context(params: GetPathwayContextInput) -> str:
    """Retrieve biological pathway membership for a gene from Reactome.

    Use this tool to understand where a target sits in cellular biology:
    which pathways it participates in, and which other genes share those pathways.
    Useful for identifying combination therapy opportunities, resistance mechanisms,
    and on/off-target pathway effects.

    Args:
        params (GetPathwayContextInput): gene_symbol, response_format.

    Returns:
        Markdown with top enriched Reactome pathways, enrichment p-values,
        pathway categories, and gene counts per pathway.
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    result = await mcp.state.reactome.get_pathway_context(symbol)
    if result is not None and not result.pathways:
        result = None
    return _fmt(result, params.response_format, f"No Reactome pathway data found for '{symbol}'.")


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_pathway_members(params: GetPathwayMembersInput) -> str:
    """Enumerate all gene members of a named Reactome pathway.

    Use this tool when you need the full gene list for a pathway family to enable
    systematic screening — e.g. "find all MAPK kinases" or "list genes in the
    PI3K/AKT pathway". Returns HGNC gene symbols that can then be passed to
    get_drug_history, get_cancer_dependency, or prioritize_target for each member.

    Accepts either a pathway display name (fuzzy search) or a Reactome stable ID
    (e.g. R-HSA-5673001) for exact lookup.

    Args:
        params (GetPathwayMembersInput): pathway_name_or_id, max_genes.

    Returns:
        Markdown list of HGNC gene symbols for all ReferenceGeneProduct participants
        in the pathway, or an error message if the pathway cannot be found.
    """
    genes = await mcp.state.reactome.get_pathway_members(
        params.pathway_name_or_id, params.max_genes
    )
    if not genes:
        return (
            f"**No pathway members found for '{params.pathway_name_or_id}'.**\n\n"
            "Try a more specific name (e.g. 'RAF/MAP kinase cascade') or use a "
            "Reactome stable ID (e.g. R-HSA-5673001)."
        )
    lines = [
        f"## Pathway members: {params.pathway_name_or_id}",
        f"**{len(genes)} gene(s)** found in Reactome.\n",
        "| # | Gene |",
        "|---|------|",
    ]
    for i, gene in enumerate(genes, 1):
        lines.append(f"| {i} | {gene} |")
    return "\n".join(lines)


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=False, openWorldHint=True
    )
)
async def prioritize_target(params: PrioritizeTargetInput) -> str:
    """Generate a comprehensive drug discovery target assessment report.

    Use this tool when you need a full evidence synthesis for target prioritization.
    Queries all databases in parallel (UniProt, Open Targets, DepMap, GWAS Catalog,
    PubChem, ChEMBL) and returns a structured Markdown report with a composite
    priority score (0–10). Use for go/no-go decisions on a target–indication pair.

    Set extended=True for a deeper report that also includes AlphaFold structure data,
    STRING interactome (selectivity risks), drug history (DGIdb + ClinicalTrials.gov),
    and Reactome pathway context. Extended mode adds ~10–20s to query time.

    Args:
        params (PrioritizeTargetInput): gene_symbol, indication, extended, response_format.

    Returns:
        Markdown report with priority_score (0–10), priority_tier (High/Medium/Low),
        evidence summary, scoring breakdown table, and data gaps.
        With extended=True, also includes structure, interactors, drugs, and pathways.
    """
    result = await _prioritize_target(
        gene_symbol=params.gene_symbol,
        indication=params.indication,
        uniprot=mcp.state.uniprot,
        open_targets=mcp.state.open_targets,
        depmap=mcp.state.depmap,
        gwas=mcp.state.gwas,
        pubchem=mcp.state.pubchem,
        chembl=mcp.state.chembl,
        alphafold=mcp.state.alphafold if params.extended else None,
        string_db=mcp.state.string_db if params.extended else None,
        dgidb=mcp.state.dgidb if params.extended else None,
        clinical_trials=mcp.state.clinical_trials if params.extended else None,
        reactome=mcp.state.reactome if params.extended else None,
    )
    return _fmt(result, params.response_format, "")


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=False, openWorldHint=True
    )
)
async def compare_targets(params: CompareTargetsInput) -> str:
    """Compare 2–5 drug targets side by side for a given therapeutic indication.

    Runs a full prioritize_target assessment for each gene in parallel and returns
    a ranked Markdown comparison table. Use this tool when you need to choose between
    multiple candidate targets or rank a target panel for an indication.

    Args:
        params (CompareTargetsInput): gene_symbols (2–5), indication, response_format.

    Returns:
        Markdown table ranking all genes by priority score, with per-gene evidence summaries.
    """
    gene_symbols = list(params.gene_symbols)
    dropped: list[str] = []
    if len(gene_symbols) > 5:
        dropped = gene_symbols[5:]
        gene_symbols = gene_symbols[:5]
        logger.warning("compare_targets: capped to 5 genes, dropped: %s", dropped)

    reports = await asyncio.gather(
        *[
            _prioritize_target(
                gene_symbol=sym,
                indication=params.indication,
                uniprot=mcp.state.uniprot,
                open_targets=mcp.state.open_targets,
                depmap=mcp.state.depmap,
                gwas=mcp.state.gwas,
                pubchem=mcp.state.pubchem,
                chembl=mcp.state.chembl,
            )
            for sym in gene_symbols
        ],
        return_exceptions=True,
    )

    rows: list[TargetComparisonRow] = []
    for sym, report in zip(gene_symbols, reports):
        if isinstance(report, Exception):
            logger.warning("compare_targets: failed for %s — %s", sym, report)
            rows.append(
                TargetComparisonRow(
                    gene_symbol=sym.upper(),
                    priority_score=0.0,
                    priority_tier="Error",
                    data_gaps=["all"],
                    evidence_summary=f"Query failed: {report}",
                )
            )
            continue

        cd = report.cancer_dependency
        depmap_pct = int(cd.fraction_dependent_lines * 100) if cd else None
        depmap_real = cd is not None and "DepMap Chronos" in cd.data_source

        rows.append(
            TargetComparisonRow(
                gene_symbol=report.gene_symbol,
                priority_score=report.priority_score,
                priority_tier=report.priority_tier,
                ot_score=report.disease_association.overall_score
                if report.disease_association
                else None,
                depmap_pct=depmap_pct,
                depmap_real_data=depmap_real,
                compound_count=report.compounds.total_active_compounds
                if report.compounds
                else None,
                gwas_count=report.gwas_evidence.total_associations
                if report.gwas_evidence
                else None,
                data_gaps=report.data_gaps,
                evidence_summary=report.evidence_summary,
            )
        )

    comparison = ComparisonReport(indication=params.indication, rows=rows)
    result = _fmt(comparison, params.response_format, "")
    if dropped:
        truncation_note = (
            f"> **Note:** Input exceeded the 5-target limit. "
            f"The following gene{'s were' if len(dropped) > 1 else ' was'} dropped: "
            f"{', '.join(dropped)}. Re-run with a smaller set to include them.\n\n"
        )
        result = truncation_note + result
    return result


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=False, openWorldHint=True
    )
)
async def run_biology_workflow(params: RunBiologyWorkflowInput) -> str:
    """Answer a complex biology or drug discovery question using AI-driven tool selection.

    Use this tool for multi-step questions that require chaining several databases in a
    sequence the caller hasn't pre-determined. An internal Claude agent reasons about
    which tools to call and in what order, executes them, and synthesizes the results.

    Compared to calling individual tools directly:
    - Handles questions involving unknown pathways, novel target lists, or open-ended
      discovery tasks where the right tool sequence is not obvious in advance.
    - Dynamically adapts — if a tool returns unexpected data (e.g. a gene resolves to
      a different symbol), the agent adjusts subsequent calls accordingly.
    - Returns a synthesized Markdown answer citing specific evidence numbers.

    Args:
        params (RunBiologyWorkflowInput): question.

    Returns:
        Synthesized Markdown answer with citations from all consulted databases.
    """
    registry = build_tool_registry(mcp.state)
    return await run_agent_loop(params.question, registry)


# ---------------------------------------------------------------------------
# MCP Resource: tool discovery
# ---------------------------------------------------------------------------


@mcp.resource("tool://registry")
async def tool_registry_resource() -> str:
    """Structured reference of all available tools for agent and client discovery.

    This resource is the MCP-native discovery mechanism: clients and orchestrators
    can read it at connection time to understand the server's capabilities without
    invoking any tool.  Each tool entry includes its description, ``use_when``
    guidance (written to be embedding-searchable at scale), required inputs, and
    category grouping.

    Differs from ``list_tools`` in that it provides richer semantic metadata
    designed for tool selection reasoning, not just name and schema.

    Returns:
        Markdown document grouped by tool category with descriptions and
        ``use_when`` fields for all 13 registered tools.
    """
    registry = build_tool_registry(mcp.state)
    return format_registry_docs(registry)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
