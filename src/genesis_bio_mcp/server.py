"""genesis-bio-mcp MCP server.

Exposes 7 tools for biomedical database queries:
  - resolve_gene          UniProt + NCBI: gene symbol → canonical IDs
  - get_protein_info      UniProt Swiss-Prot protein annotation
  - get_target_disease    Open Targets: target–disease association score
  - get_cancer_dependency DepMap: CRISPR essentiality across cancer lines
  - get_gwas_evidence     GWAS Catalog: genetic associations for a trait
  - get_compounds         PubChem: active small molecules against a target
  - prioritize_target     Orchestration: full target assessment report
"""

from __future__ import annotations

import logging
from contextlib import asynccontextmanager

import httpx
from mcp.server.fastmcp import FastMCP

from genesis_bio_mcp.clients.depmap import DepMapClient
from genesis_bio_mcp.clients.gwas import GwasClient
from genesis_bio_mcp.clients.open_targets import OpenTargetsClient
from genesis_bio_mcp.clients.pubchem import PubChemClient
from genesis_bio_mcp.clients.uniprot import UniProtClient
from genesis_bio_mcp.tools.gene_resolver import resolve_gene as _resolve_gene
from genesis_bio_mcp.tools.target_prioritization import prioritize_target as _prioritize_target

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

_HEADERS = {
    "User-Agent": "genesis-bio-mcp/0.1 (research; github.com/WSobo/genesis-bio-mcp)",
    "Accept": "application/json",
}


@asynccontextmanager
async def lifespan(server: FastMCP):
    """Manage a shared httpx.AsyncClient across all tool invocations."""
    async with httpx.AsyncClient(
        headers=_HEADERS, timeout=30.0, follow_redirects=True
    ) as client:
        server.state.uniprot = UniProtClient(client)
        server.state.open_targets = OpenTargetsClient(client)
        server.state.depmap = DepMapClient(client)
        server.state.gwas = GwasClient(client)
        server.state.pubchem = PubChemClient(client)
        yield


mcp = FastMCP("genesis-bio-mcp", lifespan=lifespan)


# ---------------------------------------------------------------------------
# Tool definitions
# ---------------------------------------------------------------------------


@mcp.tool()
async def resolve_gene(gene_name: str) -> dict:
    """Resolve a gene name or alias to canonical identifiers across databases.

    Use this tool FIRST when the input is a gene alias, synonym, or informal name.
    Provides HGNC symbol, NCBI Gene ID, and UniProt accession needed by other tools.

    Args:
        gene_name: Gene symbol, alias, or synonym. Case-insensitive.
                   Examples: 'BRAF', 'braf', 'B-RAF', 'BRAF1', 'p44erk1'

    Returns:
        GeneResolution with hgnc_symbol, ncbi_gene_id, uniprot_accession, and synonyms.
    """
    result = await _resolve_gene(gene_name, uniprot_client=mcp.state.uniprot)
    return result.model_dump()


@mcp.tool()
async def get_protein_info(gene_symbol: str) -> dict:
    """Retrieve protein-level annotation for a human gene from UniProt Swiss-Prot.

    Use this tool to understand a protein's biological function, subcellular location,
    known disease-linked variants, and available 3D structures. Best used after
    resolve_gene to ensure a canonical symbol.

    Args:
        gene_symbol: Approved HGNC gene symbol (uppercase).
                     Examples: 'BRAF', 'EGFR', 'TP53', 'PCSK9'

    Returns:
        ProteinInfo with function summary, pathways, disease associations, PDB IDs,
        known variants, and reviewed status.
    """
    result = await mcp.state.uniprot.get_protein(gene_symbol)
    if result is None:
        return {
            "error": f"No UniProt Swiss-Prot entry found for gene '{gene_symbol}' in Homo sapiens."
        }
    return result.model_dump()


@mcp.tool()
async def get_target_disease_association(gene_symbol: str, disease_name: str) -> dict:
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
        TargetDiseaseAssociation with overall_score and per-datatype evidence scores
        (genetic_association, somatic_mutation, known_drug, literature).
    """
    result = await mcp.state.open_targets.get_association(gene_symbol, disease_name)
    if result is None:
        return {
            "error": f"No Open Targets association found for '{gene_symbol}' and '{disease_name}'."
        }
    return result.model_dump()


@mcp.tool()
async def get_cancer_dependency(gene_symbol: str) -> dict:
    """Retrieve CRISPR essentiality scores for a gene across cancer cell lines from DepMap.

    Use this tool to assess whether a gene is a cancer dependency — i.e., cancer cells
    require it for survival (CERES score < -0.5). Pan-essential genes (>90% of lines
    dependent) are core cellular machinery and may have narrow therapeutic windows.

    Args:
        gene_symbol: HGNC gene symbol. Examples: 'BRAF', 'MYC', 'CDK4', 'EGFR'

    Returns:
        CancerDependency with mean CERES score, fraction of dependent lines, pan-essential
        flag, and cancer lineage breakdown.
    """
    result = await mcp.state.depmap.get_essentiality(gene_symbol)
    if result is None:
        return {
            "error": (
                f"DepMap essentiality data unavailable for '{gene_symbol}'. "
                f"Check manually at https://depmap.org/portal/gene/{gene_symbol}"
            )
        }
    return result.model_dump()


@mcp.tool()
async def get_gwas_evidence(gene_symbol: str, trait: str) -> dict:
    """Retrieve GWAS Catalog associations linking a gene locus to a disease trait.

    Use this tool to find genome-wide significant SNP associations (p < 5e-8) near
    a gene for a phenotypic trait. High association counts and low p-values strengthen
    genetic causality evidence for target selection.

    Args:
        gene_symbol: HGNC gene symbol. Examples: 'FTO', 'APOE', 'BRCA1', 'PCSK9'
        trait: Trait or disease name for filtering. Case-insensitive substring match.
               Examples: 'type 2 diabetes', 'body mass index', 'breast cancer', 'melanoma'

    Returns:
        GwasEvidence with list of GWAS hits including p-values, effect sizes, and study IDs.
    """
    result = await mcp.state.gwas.get_evidence(gene_symbol, trait)
    if result is None:
        return {
            "error": f"No GWAS Catalog associations found for gene '{gene_symbol}' and trait '{trait}'."
        }
    return result.model_dump()


@mcp.tool()
async def get_compounds(gene_symbol: str) -> dict:
    """Retrieve small molecules with bioactivity against a gene target from PubChem.

    Use this tool to assess target druggability — whether active chemical matter exists.
    Returns active compounds sorted by potency (lowest IC50/EC50 first). A count >50
    indicates a well-explored chemical space and tractable target.

    Args:
        gene_symbol: HGNC gene symbol. Examples: 'BRAF', 'EGFR', 'CDK2', 'HDAC1'

    Returns:
        Compounds with total active compound count and top 20 by potency (IC50/EC50 in nM).
    """
    result = await mcp.state.pubchem.get_compounds(gene_symbol)
    if result is None:
        return {
            "error": f"No PubChem bioactivity data found for gene '{gene_symbol}'."
        }
    return result.model_dump()


@mcp.tool()
async def prioritize_target(gene_symbol: str, indication: str) -> dict:
    """Generate a comprehensive drug discovery target assessment report.

    Use this tool when you need a full evidence synthesis for target prioritization.
    Queries all databases in parallel (UniProt, Open Targets, DepMap, GWAS Catalog,
    PubChem) and returns a structured report with a composite priority score (0–10).
    Use for go/no-go decisions on a target–indication pair.

    Args:
        gene_symbol: HGNC gene symbol of the candidate drug target. Example: 'BRAF'
        indication: Therapeutic indication or disease area.
                    Examples: 'melanoma', 'non-small cell lung cancer', 'colorectal cancer'

    Returns:
        TargetPrioritizationReport with priority_score (0–10), priority_tier
        (High/Medium/Low), evidence_summary, full evidence from all databases,
        and data_gaps listing any sources that could not be retrieved.
    """
    result = await _prioritize_target(
        gene_symbol=gene_symbol,
        indication=indication,
        uniprot=mcp.state.uniprot,
        open_targets=mcp.state.open_targets,
        depmap=mcp.state.depmap,
        gwas=mcp.state.gwas,
        pubchem=mcp.state.pubchem,
    )
    return result.model_dump()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
