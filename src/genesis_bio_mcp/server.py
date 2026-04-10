"""genesis-bio-mcp MCP server.

Exposes 12 tools for biomedical database queries:
  - resolve_gene             UniProt + NCBI: gene symbol → canonical IDs
  - get_protein_info         UniProt Swiss-Prot protein annotation
  - get_target_disease       Open Targets: target–disease association score
  - get_cancer_dependency    DepMap: CRISPR essentiality across cancer lines
  - get_gwas_evidence        GWAS Catalog: genetic associations for a trait
  - get_compounds            PubChem: active small molecules against a target
  - get_protein_structure    AlphaFold + RCSB PDB: structural data
  - get_protein_interactome  STRING: binding partners and selectivity risks
  - get_drug_history         DGIdb + ClinicalTrials.gov: known drugs and trials
  - get_pathway_context      Reactome: pathway membership and enrichment
  - prioritize_target        Orchestration: full target assessment report
  - compare_targets          Compare 2–5 targets side by side for an indication

All tools return Markdown strings for direct LLM consumption.
"""

from __future__ import annotations

import asyncio
import logging
from contextlib import asynccontextmanager

import httpx
from mcp.server.fastmcp import FastMCP

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
from genesis_bio_mcp.models import ComparisonReport, DrugHistory, TargetComparisonRow
from genesis_bio_mcp.tools.gene_resolver import resolve_gene as _resolve_gene
from genesis_bio_mcp.tools.target_prioritization import (
    prioritize_target as _prioritize_target,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

_HEADERS = {
    "User-Agent": "genesis-bio-mcp/0.1 (research; github.com/WSobo/genesis-bio-mcp)",
    "Accept": "application/json",
}


@asynccontextmanager
async def lifespan(server: FastMCP):
    """Manage a shared httpx.AsyncClient and pre-load the DepMap gene cache."""
    async with httpx.AsyncClient(headers=_HEADERS, timeout=30.0, follow_redirects=True) as client:
        # Fetch DepMap gene_dep_summary once at startup for instant lookups
        gene_dep_cache = await load_depmap_cache(client)

        server.state.uniprot = UniProtClient(client)
        server.state.open_targets = OpenTargetsClient(client)
        server.state.depmap = DepMapClient(client, gene_dep_cache)
        server.state.gwas = GwasClient(client)
        server.state.pubchem = PubChemClient(client)
        server.state.chembl = ChEMBLClient(client)
        server.state.alphafold = AlphaFoldClient(client)
        server.state.string_db = StringDbClient(client)
        server.state.dgidb = DGIdbClient(client)
        server.state.clinical_trials = ClinicalTrialsClient(client)
        server.state.reactome = ReactomeClient(client)
        yield


mcp = FastMCP("genesis-bio-mcp", lifespan=lifespan)


# ---------------------------------------------------------------------------
# Tool definitions — all return Markdown strings for LLM readability
# ---------------------------------------------------------------------------


@mcp.tool()
async def resolve_gene(gene_name: str) -> str:
    """Resolve a gene name or alias to canonical identifiers across databases.

    Use this tool FIRST when the input is a gene alias, synonym, or informal name.
    Provides HGNC symbol, NCBI Gene ID, and UniProt accession needed by other tools.

    Args:
        gene_name: Gene symbol, alias, or synonym. Case-insensitive.
                   Examples: 'BRAF', 'braf', 'B-RAF', 'BRAF1', 'p44erk1'

    Returns:
        Markdown with canonical symbol, NCBI Gene ID, UniProt accession, and synonyms.
    """
    result = await _resolve_gene(gene_name, uniprot_client=mcp.state.uniprot)
    return result.to_markdown()


@mcp.tool()
async def get_protein_info(gene_symbol: str) -> str:
    """Retrieve protein-level annotation for a human gene from UniProt Swiss-Prot.

    Use this tool to understand a protein's biological function, subcellular location,
    known disease-linked variants, and available 3D structures. Best used after
    resolve_gene to ensure a canonical symbol.

    Args:
        gene_symbol: Approved HGNC gene symbol (uppercase).
                     Examples: 'BRAF', 'EGFR', 'TP53', 'PCSK9'

    Returns:
        Markdown with function summary, pathways, disease associations, PDB IDs,
        known variants, and reviewed status.
    """
    result = await mcp.state.uniprot.get_protein(gene_symbol)
    if result is None:
        return f"**Error:** No UniProt Swiss-Prot entry found for gene '{gene_symbol}' in Homo sapiens."
    return result.to_markdown()


@mcp.tool()
async def get_target_disease_association(gene_symbol: str, disease_name: str) -> str:
    """Query Open Targets for the evidence-based association score between a gene and disease.

    Use this tool to assess genetic, clinical, and literature evidence linking a drug
    target to a specific indication. The overall_score (0–1) is Open Targets' aggregate
    evidence strength; scores >0.5 are considered strong support for a target–disease link.

    Args:
        gene_symbol: HGNC gene symbol. Example: 'BRAF'
        disease_name: Free-text disease or indication name. Open Targets maps this
                      to the closest EFO ontology term automatically.
                      Examples: 'melanoma', 'non-small cell lung cancer', 'Alzheimer disease'

    Returns:
        Markdown with overall_score and per-datatype evidence scores
        (genetic_association, somatic_mutation, known_drug, literature).
    """
    result = await mcp.state.open_targets.get_association(gene_symbol, disease_name)
    if result is None:
        return f"**Error:** No Open Targets association found for '{gene_symbol}' and '{disease_name}'."
    return result.to_markdown()


@mcp.tool()
async def get_cancer_dependency(gene_symbol: str) -> str:
    """Retrieve CRISPR essentiality scores for a gene across cancer cell lines from DepMap.

    Uses real DepMap Chronos Combined data when available (loaded at server startup),
    supplemented by Open Targets somatic mutation data for lineage context.
    Pan-essential genes (common_essential=True) are core cellular machinery and may
    have narrow therapeutic windows.

    Args:
        gene_symbol: HGNC gene symbol. Examples: 'BRAF', 'MYC', 'CDK4', 'EGFR'

    Returns:
        Markdown with fraction of dependent lines, pan-essential flag, top lineages,
        and the data source (real DepMap or OT proxy).
    """
    result = await mcp.state.depmap.get_essentiality(gene_symbol)
    if result is None:
        return (
            f"**Error:** DepMap essentiality data unavailable for '{gene_symbol}'. "
            f"Check manually at https://depmap.org/portal/gene/{gene_symbol}"
        )
    return result.to_markdown()


@mcp.tool()
async def get_gwas_evidence(gene_symbol: str, trait: str) -> str:
    """Retrieve GWAS Catalog associations linking a gene locus to a disease trait.

    Use this tool to find genome-wide significant SNP associations (p < 5e-8) near
    a gene for a phenotypic trait. High association counts and low p-values strengthen
    genetic causality evidence for target selection.

    Args:
        gene_symbol: HGNC gene symbol. Examples: 'FTO', 'APOE', 'BRCA1', 'PCSK9'
        trait: Trait or disease name for filtering. Case-insensitive substring match.
               Examples: 'type 2 diabetes', 'body mass index', 'breast cancer', 'melanoma'

    Returns:
        Markdown with GWAS hit count, strongest p-value, and top associations table.
    """
    result = await mcp.state.gwas.get_evidence(gene_symbol, trait)
    if result is None:
        return f"**Error:** No GWAS Catalog associations found for gene '{gene_symbol}' and trait '{trait}'."
    return result.to_markdown()


@mcp.tool()
async def get_compounds(gene_symbol: str) -> str:
    """Retrieve small molecules with bioactivity against a gene target from PubChem.

    Use this tool to assess target druggability — whether active chemical matter exists.
    Returns active compounds sorted by potency (lowest IC50/EC50 first). A count >50
    indicates a well-explored chemical space and tractable target.

    Args:
        gene_symbol: HGNC gene symbol. Examples: 'BRAF', 'EGFR', 'CDK2', 'HDAC1'

    Returns:
        Markdown with total active compound count and top compounds by potency (IC50/EC50 in nM).
    """
    result = await mcp.state.pubchem.get_compounds(gene_symbol)
    if result is None:
        return f"**Error:** No PubChem bioactivity data found for gene '{gene_symbol}'."
    return result.to_markdown()


@mcp.tool()
async def get_protein_structure(gene_symbol: str) -> str:
    """Retrieve structural data for a protein from AlphaFold and RCSB PDB.

    Use this tool to assess structural feasibility for drug design. Reports AlphaFold
    prediction confidence (pLDDT), experimental PDB structures, best resolution,
    and whether inhibitor-bound structures exist (critical for structure-based design).

    Args:
        gene_symbol: HGNC gene symbol. Examples: 'BRAF', 'EGFR', 'CDK2', 'KRAS'

    Returns:
        Markdown with AlphaFold pLDDT score, PDB structure count and resolution,
        and whether ligand-bound structures are available.
    """
    # Resolve UniProt accession first
    protein = await mcp.state.uniprot.get_protein(gene_symbol)
    accession = protein.uniprot_accession if protein else None
    result = await mcp.state.alphafold.get_structure(gene_symbol, uniprot_accession=accession)
    if result is None:
        return f"**Error:** No structural data found for '{gene_symbol}'. Ensure the gene symbol is correct."
    return result.to_markdown()


@mcp.tool()
async def get_protein_interactome(gene_symbol: str) -> str:
    """Retrieve protein interaction partners from STRING to assess selectivity risks.

    Use this tool to understand which proteins interact with your target — binding
    partners, paralogs, and pathway neighbors that your compound might also engage.
    High-confidence interactors (score ≥ 0.9) are particularly important for selectivity.

    Args:
        gene_symbol: HGNC gene symbol. Examples: 'BRAF', 'EGFR', 'CDK2', 'JAK2'

    Returns:
        Markdown with top 20 interaction partners sorted by STRING confidence score,
        evidence types (experiments, database, coexpression), and total partner count.
    """
    result = await mcp.state.string_db.get_interactome(gene_symbol)
    if result is None or result.total_partners == 0:
        return f"**Error:** No STRING interaction data found for '{gene_symbol}'."
    return result.to_markdown()


@mcp.tool()
async def get_drug_history(gene_symbol: str) -> str:
    """Retrieve the drug development history for a gene target.

    Use this tool to understand the clinical precedent for targeting a gene:
    what drugs already exist (approved or investigational), what indications have
    been pursued, and how many clinical trials have targeted this gene. Essential
    for first-in-class vs. best-in-class strategy decisions.

    Args:
        gene_symbol: HGNC gene symbol. Examples: 'BRAF', 'EGFR', 'PCSK9', 'KRAS'

    Returns:
        Markdown with known drugs (type, phase, approval status), trial counts by
        phase, and a table of recent clinical trials from ClinicalTrials.gov.
    """
    drugs, ct_result = await asyncio.gather(
        mcp.state.dgidb.get_drug_interactions(gene_symbol),
        mcp.state.clinical_trials.get_trials(gene_symbol),
    )
    ct_trials, ct_counts = ct_result
    approved_count = sum(1 for d in drugs if d.approved)
    result = DrugHistory(
        gene_symbol=gene_symbol,
        known_drugs=drugs,
        approved_drug_count=approved_count,
        trial_counts_by_phase=ct_counts,
        recent_trials=ct_trials[:10],
    )
    if not drugs and not ct_trials:
        return f"**No drug history found for '{gene_symbol}'.** This may be a first-in-class opportunity."
    return result.to_markdown()


@mcp.tool()
async def get_pathway_context(gene_symbol: str) -> str:
    """Retrieve biological pathway membership for a gene from Reactome.

    Use this tool to understand where a target sits in cellular biology:
    which pathways it participates in, and which other genes share those pathways.
    Useful for identifying combination therapy opportunities, resistance mechanisms,
    and on/off-target pathway effects.

    Args:
        gene_symbol: HGNC gene symbol. Examples: 'BRAF', 'EGFR', 'MTOR', 'STAT3'

    Returns:
        Markdown with top enriched Reactome pathways, enrichment p-values,
        pathway categories, and gene counts per pathway.
    """
    result = await mcp.state.reactome.get_pathway_context(gene_symbol)
    if result is None or not result.pathways:
        return f"**Error:** No Reactome pathway data found for '{gene_symbol}'."
    return result.to_markdown()


@mcp.tool()
async def prioritize_target(gene_symbol: str, indication: str, extended: bool = False) -> str:
    """Generate a comprehensive drug discovery target assessment report.

    Use this tool when you need a full evidence synthesis for target prioritization.
    Queries all databases in parallel (UniProt, Open Targets, DepMap, GWAS Catalog,
    PubChem, ChEMBL) and returns a structured Markdown report with a composite
    priority score (0–10). Use for go/no-go decisions on a target–indication pair.

    Set extended=True for a deeper report that also includes AlphaFold structure data,
    STRING interactome (selectivity risks), drug history (DGIdb + ClinicalTrials.gov),
    and Reactome pathway context. Extended mode adds ~10–20s to query time.

    Args:
        gene_symbol: HGNC gene symbol of the candidate drug target. Example: 'BRAF'
        indication: Therapeutic indication or disease area.
                    Examples: 'melanoma', 'non-small cell lung cancer', 'colorectal cancer'
        extended: If True, also retrieves structural, interactome, drug history, and
                  pathway data. Default False for faster standard reports.

    Returns:
        Markdown report with priority_score (0–10), priority_tier (High/Medium/Low),
        evidence summary, scoring breakdown table, and data gaps.
        With extended=True, also includes structure, interactors, drugs, and pathways.
    """
    result = await _prioritize_target(
        gene_symbol=gene_symbol,
        indication=indication,
        uniprot=mcp.state.uniprot,
        open_targets=mcp.state.open_targets,
        depmap=mcp.state.depmap,
        gwas=mcp.state.gwas,
        pubchem=mcp.state.pubchem,
        chembl=mcp.state.chembl,
        alphafold=mcp.state.alphafold if extended else None,
        string_db=mcp.state.string_db if extended else None,
        dgidb=mcp.state.dgidb if extended else None,
        clinical_trials=mcp.state.clinical_trials if extended else None,
        reactome=mcp.state.reactome if extended else None,
    )
    return result.to_markdown()


@mcp.tool()
async def compare_targets(gene_symbols: list[str], indication: str) -> str:
    """Compare 2–5 drug targets side by side for a given therapeutic indication.

    Runs a full prioritize_target assessment for each gene in parallel and returns
    a ranked Markdown comparison table. Use this tool when you need to choose between
    multiple candidate targets or rank a target panel for an indication.

    Args:
        gene_symbols: List of 2–5 HGNC gene symbols. Example: ['BRAF', 'EGFR', 'KRAS']
        indication: Therapeutic indication shared across all targets.
                    Example: 'melanoma', 'non-small cell lung cancer'

    Returns:
        Markdown table ranking all genes by priority score, with per-gene evidence summaries.
    """
    if len(gene_symbols) < 2:
        return "**Error:** compare_targets requires at least 2 gene symbols."
    if len(gene_symbols) > 5:
        gene_symbols = gene_symbols[:5]
        logger.warning("compare_targets: capped to 5 genes")

    reports = await asyncio.gather(
        *[
            _prioritize_target(
                gene_symbol=sym,
                indication=indication,
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

    comparison = ComparisonReport(indication=indication, rows=rows)
    return comparison.to_markdown()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
