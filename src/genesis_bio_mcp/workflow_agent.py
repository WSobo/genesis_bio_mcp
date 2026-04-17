"""AI-powered workflow agent for dynamic tool selection and multi-step chaining.

Provides:
  - ToolSpec: typed description of a callable tool for agent reasoning
  - build_tool_registry: constructs a callable registry from live server state
  - run_agent_loop: Claude-powered agentic loop that dynamically selects and
    chains tools to answer free-text biology questions
"""

from __future__ import annotations

import asyncio
import logging
from collections.abc import Awaitable, Callable
from dataclasses import dataclass
from typing import Any

import anthropic

from genesis_bio_mcp.config.settings import settings
from genesis_bio_mcp.models import ComparisonReport, DrugHistory, TargetComparisonRow
from genesis_bio_mcp.tools.gene_resolver import resolve_gene as _resolve_gene
from genesis_bio_mcp.tools.target_prioritization import prioritize_target as _prioritize_target

logger = logging.getLogger(__name__)

_MODEL = settings.claude_model
_MAX_TOKENS = 8192
_SYSTEM_PREAMBLE = """\
You are a drug discovery research agent with access to biomedical databases.
Answer the user's biology question by calling the available tools in a logical order.

Workflow guidance:
1. Start with resolve_gene if the gene name might be an alias or informal name.
2. Use pathway/annotation tools (get_pathway_context, get_protein_info) for broad biology questions. Use get_pathway_members to enumerate all genes in a pathway for systematic screening.
3. Use evidence tools (get_target_disease_association, get_cancer_dependency, get_gwas_evidence) for disease links.
4. Use druggability tools (get_compounds, get_drug_history) to assess chemical tractability and competition.
5. Call prioritize_target for a full evidence synthesis and composite score on a specific target–indication pair.
6. Call compare_targets when ranking 2–5 candidates side by side.
7. For antibody/nanobody design questions: call get_antibody_structures to find PDB-curated structural templates; use get_protein_structure for the antigen conformation.
8. For protein engineering questions: call get_variant_constraints first to understand mutation tolerance (pLI/LOEUF); then use get_protein_structure for structural context.

Always synthesize the tool results into a coherent Markdown answer. \
Cite specific numbers (scores, counts, p-values) from tool outputs.\
"""


@dataclass
class ToolSpec:
    """Typed description of a callable tool for agent-driven selection.

    Fields are designed for two complementary retrieval modes:

    1. Coarse filter — ``tool_category`` partitions tools into a small set of
       buckets (gene_annotation, disease_evidence, druggability, structure,
       pathways, synthesis).  A category filter cuts candidate tools from N to
       ~2–4 before any LLM reasoning is needed.

    2. Semantic retrieval — ``use_when`` is a free-text sentence written
       specifically to be **embedding-searchable**.  At registry scale (50+
       tools), embed all ``use_when`` strings once at startup and retrieve
       relevant tools per query via cosine similarity:

       .. code-block:: python

           def find_tools(
               query: str,
               registry: dict[str, "ToolSpec"],
               *,
               embedder: Callable[[str], list[float]],
               top_k: int = 4,
           ) -> list["ToolSpec"]:
               q_vec = embedder(query)
               scores = {
                   name: cosine(q_vec, embedder(spec.use_when))
                   for name, spec in registry.items()
               }
               return [registry[n] for n in sorted(scores, key=scores.__getitem__, reverse=True)[:top_k]]

       Passing only the top-k specs to Claude instead of all N reduces prompt
       size, lowers latency, and decreases hallucinated tool names on large
       registries.
    """

    name: str
    description: str
    input_schema: dict  # JSON Schema for Claude tool_use input
    tool_category: str  # coarse pre-filter bucket before semantic retrieval
    use_when: str  # embedding-searchable sentence describing when to invoke this tool
    fn: Callable[..., Awaitable[str]]


def build_tool_registry(state: Any) -> dict[str, ToolSpec]:
    """Construct a live tool registry from initialized server state.

    Each ToolSpec.fn wraps an existing client method so the agent can call
    them directly by name without going through the MCP protocol layer.

    Args:
        state: FastMCP server state object with all database clients initialized.

    Returns:
        Mapping of tool name → ToolSpec with a callable async fn.
    """
    # ---- async wrappers for tools that need special handling ----

    async def _resolve_gene_fn(gene_name: str) -> str:
        result = await _resolve_gene(gene_name, uniprot_client=state.uniprot)
        return result.to_markdown()

    async def _get_protein_info_fn(gene_symbol: str) -> str:
        result = await state.uniprot.get_protein(gene_symbol)
        if result is None:
            return f"No UniProt Swiss-Prot entry found for '{gene_symbol}'."
        return result.to_markdown()

    async def _get_protein_sequence_fn(
        gene_symbol: str, start: int | None = None, end: int | None = None
    ) -> str:
        from genesis_bio_mcp.models import ProteinSequence
        from genesis_bio_mcp.tools.biochem import compute_features, scan_liabilities

        protein = await state.uniprot.get_protein(gene_symbol)
        if protein is None:
            return f"No UniProt Swiss-Prot entry found for '{gene_symbol}'."
        fetched = await state.uniprot.get_sequence(protein.uniprot_accession, start=start, end=end)
        if fetched is None:
            return (
                f"Could not fetch FASTA for '{gene_symbol}' (UniProt {protein.uniprot_accession})."
            )
        sequence, organism, description = fetched
        if start is not None and end is not None:
            offset = start - 1
            local_disulfide = {
                p - offset for p in protein.disulfide_bond_positions if start <= p <= end
            }
            region_disulfide = sorted(local_disulfide)
        else:
            local_disulfide = set(protein.disulfide_bond_positions)
            region_disulfide = list(protein.disulfide_bond_positions)
        result = ProteinSequence(
            uniprot_accession=protein.uniprot_accession,
            gene_symbol=protein.gene_symbol,
            organism=organism or protein.organism,
            description=description or protein.protein_name,
            sequence=sequence,
            region_start=start,
            region_end=end,
            features=compute_features(sequence),
            liabilities=scan_liabilities(
                sequence, disulfide_annotated_positions=local_disulfide or None
            ),
            disulfide_bond_positions=region_disulfide,
        )
        return result.to_markdown()

    async def _get_target_disease_fn(gene_symbol: str, disease_name: str) -> str:
        result = await state.open_targets.get_association(gene_symbol, disease_name)
        if result is None:
            return f"No Open Targets association found for '{gene_symbol}' and '{disease_name}'."
        return result.to_markdown()

    async def _get_cancer_dependency_fn(gene_symbol: str) -> str:
        result = await state.depmap.get_essentiality(gene_symbol)
        if result is None:
            return f"No DepMap essentiality data found for '{gene_symbol}'."
        return result.to_markdown()

    async def _get_gwas_evidence_fn(gene_symbol: str, trait: str) -> str:
        result = await state.gwas.get_evidence(gene_symbol, trait)
        if result is None:
            return f"No GWAS Catalog associations found for '{gene_symbol}' and '{trait}'."
        return result.to_markdown()

    async def _get_compounds_fn(gene_symbol: str) -> str:
        result = await state.pubchem.get_compounds(gene_symbol)
        if result is None:
            return f"No PubChem bioactivity data found for '{gene_symbol}'."
        return result.to_markdown()

    async def _get_protein_structure_fn(gene_symbol: str) -> str:
        protein = await state.uniprot.get_protein(gene_symbol)
        accession = protein.uniprot_accession if protein else None
        result = await state.alphafold.get_structure(gene_symbol, uniprot_accession=accession)
        if result is None:
            return f"No structural data found for '{gene_symbol}'."
        return result.to_markdown()

    async def _get_protein_interactome_fn(gene_symbol: str) -> str:
        result = await state.string_db.get_interactome(gene_symbol)
        if result is None or result.total_partners == 0:
            return f"No STRING interaction data found for '{gene_symbol}'."
        return result.to_markdown()

    async def _get_biogrid_interactions_fn(gene_symbol: str) -> str:
        import os

        if not os.environ.get("BIOGRID_ACCESS_KEY"):
            return (
                "BioGRID data unavailable: BIOGRID_ACCESS_KEY is not set. "
                "Register at https://webservice.thebiogrid.org/ and set the env var."
            )
        result = await state.biogrid.get_interactions(gene_symbol)
        if result is None:
            return f"No BioGRID interaction data found for '{gene_symbol}'."
        return result.to_markdown()

    async def _get_epitope_data_fn(antigen_query: str) -> str:
        result = await state.iedb.get_epitopes(antigen_query)
        if result is None:
            return f"IEDB data temporarily unavailable for '{antigen_query}'."
        if result.total_assays == 0:
            return f"No B-cell epitope data found in IEDB for '{antigen_query}'."
        return result.to_markdown()

    async def _get_antibody_structures_fn(antigen_query: str, max_results: int = 20) -> str:
        result = await state.sabdab.get_antibody_structures(antigen_query, max_results=max_results)
        if result is None:
            return f"SAbDab data temporarily unavailable for '{antigen_query}'."
        if result.total_structures == 0:
            return f"No antibody or nanobody structures found in SAbDab for '{antigen_query}'."
        return result.to_markdown()

    async def _get_variant_constraints_fn(gene_symbol: str) -> str:
        result = await state.gnomad.get_constraint(gene_symbol)
        if result is None:
            return f"No gnomAD constraint data found for '{gene_symbol}'."
        return result.to_markdown()

    async def _get_variant_effects_fn(gene_symbol: str, mutation: str) -> str:
        try:
            result = await state.variant_effects.get_effects(gene_symbol, mutation)
        except ValueError as exc:
            return f"Could not parse mutation '{mutation}': {exc}"
        return result.to_markdown()

    async def _get_domain_annotation_fn(gene_symbol: str) -> str:
        protein = await state.uniprot.get_protein(gene_symbol)
        accession = protein.uniprot_accession if protein else gene_symbol
        result = await state.interpro.get_domains(gene_symbol, accession)
        if result is None:
            return f"No InterPro domain data found for '{gene_symbol}'."
        return result.to_markdown()

    async def _get_dms_scores_fn(gene_symbol: str) -> str:
        result = await state.mavedb.get_dms_scores(gene_symbol)
        if result is None:
            return f"MaveDB data temporarily unavailable for '{gene_symbol}'."
        if result.total_score_sets == 0:
            return (
                f"No DMS score sets found in MaveDB for '{gene_symbol}'. "
                "DMS data is sparse — not all genes have been profiled."
            )
        return result.to_markdown()

    async def _get_drug_history_fn(gene_symbol: str) -> str:
        drugs, ct_result = await asyncio.gather(
            state.dgidb.get_drug_interactions(gene_symbol),
            state.clinical_trials.get_trials(gene_symbol),
        )
        ct_trials, ct_counts = ct_result
        result = DrugHistory(
            gene_symbol=gene_symbol,
            known_drugs=drugs,
            approved_drug_count=sum(1 for d in drugs if d.approved),
            trial_counts_by_phase=ct_counts,
            recent_trials=ct_trials[:10],
        )
        if not drugs and not ct_trials:
            return f"No drug history found for '{gene_symbol}'. This may be a first-in-class opportunity."
        return result.to_markdown()

    async def _get_chembl_compounds_fn(gene_symbol: str) -> str:
        result = await state.chembl.get_compounds(gene_symbol)
        if result is None:
            return f"No ChEMBL bioactivity data found for '{gene_symbol}'."
        return result.to_markdown()

    async def _get_pathway_context_fn(gene_symbol: str) -> str:
        result = await state.reactome.get_pathway_context(gene_symbol)
        if result is None or not result.pathways:
            return f"No Reactome pathway data found for '{gene_symbol}'."
        return result.to_markdown()

    async def _get_pathway_members_fn(pathway_name_or_id: str, max_genes: int = 50) -> str:
        genes = await state.reactome.get_pathway_members(pathway_name_or_id, max_genes)
        if not genes:
            return (
                f"No pathway members found for '{pathway_name_or_id}'. "
                "Try a more specific name or a Reactome stable ID (e.g. R-HSA-5673001)."
            )
        return f"Pathway members for '{pathway_name_or_id}': {', '.join(genes)}"

    async def _prioritize_target_fn(
        gene_symbol: str, indication: str, extended: bool = False
    ) -> str:
        result = await _prioritize_target(
            gene_symbol=gene_symbol,
            indication=indication,
            uniprot=state.uniprot,
            open_targets=state.open_targets,
            depmap=state.depmap,
            gwas=state.gwas,
            pubchem=state.pubchem,
            chembl=state.chembl,
            alphafold=state.alphafold if extended else None,
            string_db=state.string_db if extended else None,
            dgidb=state.dgidb if extended else None,
            clinical_trials=state.clinical_trials if extended else None,
            reactome=state.reactome if extended else None,
        )
        return result.to_markdown()

    async def _compare_targets_fn(gene_symbols: list[str], indication: str) -> str:
        if len(gene_symbols) < 2:
            return "compare_targets requires at least 2 gene symbols."
        symbols = gene_symbols[:5]
        reports = await asyncio.gather(
            *[
                _prioritize_target(
                    gene_symbol=sym,
                    indication=indication,
                    uniprot=state.uniprot,
                    open_targets=state.open_targets,
                    depmap=state.depmap,
                    gwas=state.gwas,
                    pubchem=state.pubchem,
                    chembl=state.chembl,
                )
                for sym in symbols
            ],
            return_exceptions=True,
        )
        rows: list[TargetComparisonRow] = []
        for sym, report in zip(symbols, reports):
            if isinstance(report, Exception):
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
            rows.append(
                TargetComparisonRow(
                    gene_symbol=report.gene_symbol,
                    priority_score=report.priority_score,
                    priority_tier=report.priority_tier,
                    ot_score=report.disease_association.overall_score
                    if report.disease_association
                    else None,
                    depmap_pct=int(cd.fraction_dependent_lines * 100) if cd else None,
                    depmap_real_data=cd is not None and "DepMap Chronos" in cd.data_source,
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
        return ComparisonReport(indication=indication, rows=rows).to_markdown()

    # ---- tool registry ----
    return {
        "resolve_gene": ToolSpec(
            name="resolve_gene",
            description=(
                "Resolve a gene name or alias to the canonical HGNC symbol and database IDs "
                "(NCBI Gene ID, UniProt accession). Handles aliases like HER2→ERBB2, p53→TP53."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "gene_name": {
                        "type": "string",
                        "description": "Gene symbol, alias, or synonym. Case-insensitive. Example: 'HER2', 'BRAF', 'p44erk1'",
                    }
                },
                "required": ["gene_name"],
            },
            tool_category="gene_annotation",
            use_when="Use when the input gene name might be an alias or informal name, or before any other tool to confirm the canonical symbol.",
            fn=_resolve_gene_fn,
        ),
        "get_protein_info": ToolSpec(
            name="get_protein_info",
            description=(
                "Retrieve protein annotation from UniProt Swiss-Prot: function, subcellular location, "
                "known disease variants, and PDB structures. Reviewed (Swiss-Prot) entries only."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "gene_symbol": {
                        "type": "string",
                        "description": "HGNC gene symbol. Example: 'BRAF'",
                    }
                },
                "required": ["gene_symbol"],
            },
            tool_category="gene_annotation",
            use_when="Use to understand a protein's function, biology, and disease relevance from curated annotation.",
            fn=_get_protein_info_fn,
        ),
        "get_protein_sequence": ToolSpec(
            name="get_protein_sequence",
            description=(
                "Fetch the protein FASTA sequence from UniProt and compute biochemistry "
                "(MW, theoretical pI, GRAVY, net charge at pH 7.4, ε₂₈₀ reduced/oxidized) "
                "plus a liability-motif scan (deamidation NG/NS, isomerization DG/DS, "
                "N-glycosylation NXS/NXT sequons, oxidation-prone M/W, free cysteines). "
                "Optional start/end residue range selects a sub-region. UniProt DISULFID "
                "annotations are used to distinguish bonded from free cysteines."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "gene_symbol": {
                        "type": "string",
                        "description": "HGNC gene symbol. Example: 'BRAF'",
                    },
                    "start": {
                        "type": "integer",
                        "description": "Optional 1-indexed inclusive slice start (use with end).",
                    },
                    "end": {
                        "type": "integer",
                        "description": "Optional 1-indexed inclusive slice end; must be ≥ start.",
                    },
                },
                "required": ["gene_symbol"],
            },
            tool_category="protein_engineering",
            use_when="Use to retrieve protein sequence, biochemistry, and developability liabilities before any sequence-level engineering, mutagenesis planning, or antibody CDR analysis.",
            fn=_get_protein_sequence_fn,
        ),
        "get_target_disease_association": ToolSpec(
            name="get_target_disease_association",
            description=(
                "Query Open Targets for evidence-based association score (0–1) between a gene and disease. "
                "Returns genetic_association, somatic_mutation, known_drug, and literature evidence scores."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "gene_symbol": {
                        "type": "string",
                        "description": "HGNC gene symbol. Example: 'BRAF'",
                    },
                    "disease_name": {
                        "type": "string",
                        "description": "Disease or indication. Example: 'melanoma', 'Alzheimer disease'",
                    },
                },
                "required": ["gene_symbol", "disease_name"],
            },
            tool_category="disease_evidence",
            use_when="Use to quantify the evidence strength linking a target to a specific disease. Scores >0.5 indicate strong support.",
            fn=_get_target_disease_fn,
        ),
        "get_cancer_dependency": ToolSpec(
            name="get_cancer_dependency",
            description=(
                "Retrieve CRISPR essentiality scores from DepMap: fraction of dependent cancer cell lines, "
                "pan-essential flag, and top dependent lineages. Real Chronos data preferred over OT proxy."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "gene_symbol": {
                        "type": "string",
                        "description": "HGNC gene symbol. Example: 'BRAF'",
                    }
                },
                "required": ["gene_symbol"],
            },
            tool_category="disease_evidence",
            use_when="Use to assess whether a gene is essential in cancer cell lines — key evidence for oncology targets.",
            fn=_get_cancer_dependency_fn,
        ),
        "get_gwas_evidence": ToolSpec(
            name="get_gwas_evidence",
            description=(
                "Retrieve GWAS Catalog SNP associations (p < 5e-8) near a gene for a trait. "
                "High hit counts and low p-values strengthen genetic causality evidence."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "gene_symbol": {
                        "type": "string",
                        "description": "HGNC gene symbol. Example: 'FTO'",
                    },
                    "trait": {
                        "type": "string",
                        "description": "Trait or disease name for filtering. Example: 'type 2 diabetes', 'body mass index'",
                    },
                },
                "required": ["gene_symbol", "trait"],
            },
            tool_category="disease_evidence",
            use_when="Use to find genetic evidence linking a gene locus to a disease trait via GWAS.",
            fn=_get_gwas_evidence_fn,
        ),
        "get_compounds": ToolSpec(
            name="get_compounds",
            description=(
                "Retrieve active small molecules from PubChem with bioactivity against a gene target. "
                "Returns compound count and top compounds sorted by potency (lowest IC50 first)."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "gene_symbol": {
                        "type": "string",
                        "description": "HGNC gene symbol. Example: 'BRAF'",
                    }
                },
                "required": ["gene_symbol"],
            },
            tool_category="druggability",
            use_when="Use to assess whether active chemical matter exists for a target. Count >50 indicates well-explored chemical space.",
            fn=_get_compounds_fn,
        ),
        "get_chembl_compounds": ToolSpec(
            name="get_chembl_compounds",
            description=(
                "Retrieve quantitative potency data (IC50/Ki/Kd) from ChEMBL for a gene target. "
                "Returns assay counts, pChEMBL values, and top active compounds with binding data."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "gene_symbol": {
                        "type": "string",
                        "description": "HGNC gene symbol. Example: 'BRAF'",
                    }
                },
                "required": ["gene_symbol"],
            },
            tool_category="druggability",
            use_when="Use for deeper chemical biology detail beyond PubChem: pChEMBL potency values, assay types, and binding selectivity data.",
            fn=_get_chembl_compounds_fn,
        ),
        "get_protein_structure": ToolSpec(
            name="get_protein_structure",
            description=(
                "Retrieve AlphaFold structure confidence (pLDDT) and RCSB PDB experimental structures. "
                "Reports best resolution and whether inhibitor-bound structures are available."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "gene_symbol": {
                        "type": "string",
                        "description": "HGNC gene symbol. Example: 'BRAF'",
                    }
                },
                "required": ["gene_symbol"],
            },
            tool_category="structure",
            use_when="Use to assess structural feasibility for drug design. Ligand-bound structures are critical for SBDD.",
            fn=_get_protein_structure_fn,
        ),
        "get_protein_interactome": ToolSpec(
            name="get_protein_interactome",
            description=(
                "Retrieve top STRING protein interaction partners with confidence scores. "
                "Identifies binding partners, paralogs, and pathway neighbors relevant to selectivity."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "gene_symbol": {
                        "type": "string",
                        "description": "HGNC gene symbol. Example: 'BRAF'",
                    }
                },
                "required": ["gene_symbol"],
            },
            tool_category="structure",
            use_when="Use to identify selectivity risks — proteins closely related to the target that a compound might also bind.",
            fn=_get_protein_interactome_fn,
        ),
        "get_biogrid_interactions": ToolSpec(
            name="get_biogrid_interactions",
            description=(
                "Retrieve curated protein–protein interactions from BioGRID with experimental method "
                "metadata (two-hybrid, co-IP, proximity ligation) and PubMed citations. "
                "Requires BIOGRID_ACCESS_KEY environment variable."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "gene_symbol": {
                        "type": "string",
                        "description": "HGNC gene symbol. Example: 'BRAF'",
                    }
                },
                "required": ["gene_symbol"],
            },
            tool_category="structure",
            use_when=(
                "Use for literature-curated PPI data with individual experiment records. "
                "Complements STRING when you need experimental method details or citation evidence."
            ),
            fn=_get_biogrid_interactions_fn,
        ),
        "get_epitope_data": ToolSpec(
            name="get_epitope_data",
            description=(
                "Retrieve known B-cell epitope records from IEDB for an antigen: epitope sequences, "
                "antibody isotypes, residue positions, PDB structural evidence, and publication citations. "
                "Shows which regions of the antigen surface have been characterized immunologically."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "antigen_query": {
                        "type": "string",
                        "description": (
                            "Full antigen protein name for best results. "
                            "Examples: 'epidermal growth factor receptor', 'tumor necrosis factor', "
                            "'programmed death-ligand 1'."
                        ),
                    }
                },
                "required": ["antigen_query"],
            },
            tool_category="antibody",
            use_when=(
                "Use when designing antibody therapeutics to understand which epitopes are known on the "
                "antigen, how well-characterized the target is immunologically, and whether structural "
                "epitope data (PDB) exists to guide CDR engineering."
            ),
            fn=_get_epitope_data_fn,
        ),
        "get_antibody_structures": ToolSpec(
            name="get_antibody_structures",
            description=(
                "Search SAbDab for PDB-curated antibody and nanobody (VHH) structures against a given antigen. "
                "Returns resolution, experimental method, antibody type (Fab/IgG/VHH), germline subclass, "
                "engineered flag, and binding affinity where available. Covers both conventional antibodies "
                "and camelid-derived VHH nanobodies."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "antigen_query": {
                        "type": "string",
                        "description": (
                            "Antigen gene symbol or protein name. "
                            "Examples: 'EGFR', 'HER2', 'TNF', 'programmed death-ligand 1', 'CD20'."
                        ),
                    },
                    "max_results": {
                        "type": "integer",
                        "description": "Maximum structures to return (default 20).",
                        "default": 20,
                    },
                },
                "required": ["antigen_query"],
            },
            tool_category="antibody",
            use_when=(
                "Use when designing antibody or nanobody therapeutics, finding structural templates for "
                "CDR engineering, or assessing how well-characterized an antigen is as an antibody target."
            ),
            fn=_get_antibody_structures_fn,
        ),
        "get_domain_annotation": ToolSpec(
            name="get_domain_annotation",
            description=(
                "Retrieve InterPro domain annotations for a protein: domain names, residue positions, "
                "Pfam/SMART/CDD cross-references, and GO terms. Shows which domains the protein contains "
                "and their exact boundaries."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "gene_symbol": {
                        "type": "string",
                        "description": "HGNC gene symbol. Example: 'BRAF'",
                    }
                },
                "required": ["gene_symbol"],
            },
            tool_category="protein_engineering",
            use_when=(
                "Use when planning protein engineering to understand domain boundaries — "
                "avoid mutagenesis within conserved domain cores. Complements get_variant_constraints."
            ),
            fn=_get_domain_annotation_fn,
        ),
        "get_variant_constraints": ToolSpec(
            name="get_variant_constraints",
            description=(
                "Retrieve gene-level evolutionary constraint metrics from gnomAD v4: pLI, LOEUF, oe_lof, "
                "oe_mis, Z-scores. Quantifies how much LoF and missense mutation the gene tolerates in the "
                "human population. Essential pre-filter before any protein engineering campaign."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "gene_symbol": {
                        "type": "string",
                        "description": "HGNC gene symbol. Example: 'BRAF'",
                    }
                },
                "required": ["gene_symbol"],
            },
            tool_category="protein_engineering",
            use_when=(
                "Use before designing mutations or engineering variants to understand which residues "
                "the gene tolerates losing. High pLI or low LOEUF means avoid broad mutagenesis."
            ),
            fn=_get_variant_constraints_fn,
        ),
        "get_variant_effects": ToolSpec(
            name="get_variant_effects",
            description=(
                "Retrieve mutation-level pathogenicity and fitness for a specific missense "
                "variant. Fans out to gnomAD (allele frequency, variant-ID resolution), "
                "MyVariant.info (ClinVar submissions, AlphaMissense, REVEL, CADD, SIFT, "
                "PolyPhen-2), and MaveDB (per-variant DMS fitness scores) in parallel."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "gene_symbol": {
                        "type": "string",
                        "description": "HGNC gene symbol. Example: 'TP53'",
                    },
                    "mutation": {
                        "type": "string",
                        "description": (
                            "Protein change. Accepted forms: 'R175H', 'p.R175H', "
                            "'Arg175His', 'p.Arg175His'. Example: 'R175H'"
                        ),
                    },
                },
                "required": ["gene_symbol", "mutation"],
            },
            tool_category="protein_engineering",
            use_when=(
                "Use to get the consolidated clinical and predicted-pathogenicity profile for a "
                "specific mutation. Call when the user asks 'is X variant pathogenic' or wants "
                "to check a candidate engineering mutation against ClinVar, AlphaMissense, and DMS "
                "in one step."
            ),
            fn=_get_variant_effects_fn,
        ),
        "get_dms_scores": ToolSpec(
            name="get_dms_scores",
            description=(
                "Search MaveDB for deep mutational scanning (DMS) score sets for a gene. "
                "Returns available experiments with variant counts, URNs, UniProt accessions, and citations. "
                "DMS datasets provide residue-level fitness scores for every possible amino acid change."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "gene_symbol": {
                        "type": "string",
                        "description": "HGNC gene symbol. Example: 'BRCA1'",
                    }
                },
                "required": ["gene_symbol"],
            },
            tool_category="protein_engineering",
            use_when=(
                "Use when planning specific mutations to check if experimental fitness scores exist. "
                "DMS data is the highest-fidelity residue-level tolerance signal when available. "
                "Returns empty result (not an error) when no DMS data exists for the gene."
            ),
            fn=_get_dms_scores_fn,
        ),
        "get_drug_history": ToolSpec(
            name="get_drug_history",
            description=(
                "Retrieve drug development history: known drugs from DGIdb (approved and investigational), "
                "trial counts by phase, and recent ClinicalTrials.gov entries. Essential for competitive landscape."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "gene_symbol": {
                        "type": "string",
                        "description": "HGNC gene symbol. Example: 'EGFR'",
                    }
                },
                "required": ["gene_symbol"],
            },
            tool_category="druggability",
            use_when="Use to determine clinical precedent and competitive landscape. No drugs found suggests first-in-class opportunity.",
            fn=_get_drug_history_fn,
        ),
        "get_pathway_context": ToolSpec(
            name="get_pathway_context",
            description=(
                "Retrieve Reactome pathway membership for a gene: enriched pathways with p-values, "
                "pathway categories, and participating gene counts. Useful for combination therapy reasoning."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "gene_symbol": {
                        "type": "string",
                        "description": "HGNC gene symbol. Example: 'BRAF'",
                    }
                },
                "required": ["gene_symbol"],
            },
            tool_category="pathways",
            use_when="Use to understand pathway biology, resistance mechanisms, or to find other genes in the same pathway.",
            fn=_get_pathway_context_fn,
        ),
        "get_pathway_members": ToolSpec(
            name="get_pathway_members",
            description=(
                "Enumerate all HGNC gene symbols that are members of a named Reactome pathway. "
                "Accepts a pathway display name or stable ID (R-HSA-XXXXXXX). "
                "Returns a list of gene symbols for systematic screening across a gene family."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "pathway_name_or_id": {
                        "type": "string",
                        "description": (
                            "Reactome pathway name (e.g. 'MAPK signaling', 'RAF/MAP kinase cascade') "
                            "or stable ID (e.g. 'R-HSA-5673001')."
                        ),
                    },
                    "max_genes": {
                        "type": "integer",
                        "description": "Maximum genes to return (default 50).",
                        "default": 50,
                    },
                },
                "required": ["pathway_name_or_id"],
            },
            tool_category="pathways",
            use_when=(
                "Use to enumerate all genes in a named pathway for systematic screening "
                "across a gene family, e.g. 'find all MAPK kinases' or 'list PI3K/AKT members'."
            ),
            fn=_get_pathway_members_fn,
        ),
        "prioritize_target": ToolSpec(
            name="prioritize_target",
            description=(
                "Full target prioritization report: queries all databases in parallel and returns a "
                "composite priority score (0–10) with tier (High/Medium/Low). Use for go/no-go decisions."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "gene_symbol": {
                        "type": "string",
                        "description": "HGNC gene symbol. Example: 'BRAF'",
                    },
                    "indication": {
                        "type": "string",
                        "description": "Therapeutic indication. Example: 'melanoma', 'non-small cell lung cancer'",
                    },
                    "extended": {
                        "type": "boolean",
                        "description": "If true, also retrieves structure, interactome, drug history, and pathway data. Default false.",
                    },
                },
                "required": ["gene_symbol", "indication"],
            },
            tool_category="synthesis",
            use_when="Use after gathering targeted evidence to get a full structured assessment, or as a single-call shortcut for comprehensive target analysis.",
            fn=_prioritize_target_fn,
        ),
        "compare_targets": ToolSpec(
            name="compare_targets",
            description=(
                "Compare 2–5 gene targets side by side for a given indication. "
                "Returns a ranked Markdown table with priority scores and per-gene summaries."
            ),
            input_schema={
                "type": "object",
                "properties": {
                    "gene_symbols": {
                        "type": "array",
                        "items": {"type": "string"},
                        "minItems": 2,
                        "maxItems": 5,
                        "description": "List of 2–5 HGNC gene symbols. Example: ['BRAF', 'EGFR', 'KRAS']",
                    },
                    "indication": {
                        "type": "string",
                        "description": "Shared therapeutic indication. Example: 'melanoma'",
                    },
                },
                "required": ["gene_symbols", "indication"],
            },
            tool_category="synthesis",
            use_when="Use when the question asks to rank, compare, or select among multiple candidate targets.",
            fn=_compare_targets_fn,
        ),
    }


async def run_agent_loop(
    question: str,
    registry: dict[str, ToolSpec],
    *,
    max_iterations: int = 8,
) -> str:
    """Run a Claude-powered agentic loop to answer a biology question.

    Dynamically selects and chains tools from the registry based on Claude's
    reasoning about the question. Parallel tool calls within a single turn are
    executed concurrently via asyncio.gather.

    Args:
        question: Free-text biology or drug discovery question.
        registry: Tool registry from build_tool_registry().
        max_iterations: Safety cap on agent turns (default 8).

    Returns:
        Synthesized Markdown answer.
    """
    try:
        client = anthropic.AsyncAnthropic()
    except Exception as exc:
        return (
            "**run_biology_workflow is unavailable:** `ANTHROPIC_API_KEY` is not set in the "
            "MCP server environment. Add it to `claude_desktop_config.json` under `env`, "
            "then restart Claude Desktop. For Claude Desktop users, calling individual tools "
            f"directly (get_drug_history, prioritize_target, etc.) provides equivalent capability.\n\nDetail: {exc}"
        )

    tools_for_claude = [
        {
            "name": spec.name,
            "description": (
                f"{spec.description}\n\nCategory: {spec.tool_category}. Use when: {spec.use_when}"
            ),
            "input_schema": spec.input_schema,
        }
        for spec in registry.values()
    ]

    messages: list[dict] = [{"role": "user", "content": question}]

    for iteration in range(max_iterations):
        try:
            response = await client.messages.create(
                model=_MODEL,
                max_tokens=_MAX_TOKENS,
                system=_SYSTEM_PREAMBLE,
                tools=tools_for_claude,
                messages=messages,
            )
        except anthropic.AuthenticationError:
            return (
                "**run_biology_workflow is unavailable:** `ANTHROPIC_API_KEY` is not set in the "
                "MCP server environment. Add it to `claude_desktop_config.json` under `env`, "
                "then restart Claude Desktop. For Claude Desktop users, calling individual tools "
                "directly (get_drug_history, prioritize_target, etc.) provides equivalent capability."
            )
        except Exception as exc:
            return f"**run_biology_workflow failed:** {exc}"

        if response.stop_reason == "end_turn":
            return _extract_text(response)

        if response.stop_reason == "tool_use":
            tool_use_blocks = [b for b in response.content if b.type == "tool_use"]

            # Execute all requested tools concurrently
            results = await asyncio.gather(
                *[_execute_tool(block, registry) for block in tool_use_blocks],
                return_exceptions=True,
            )

            tool_results = []
            for block, result in zip(tool_use_blocks, results):
                if isinstance(result, Exception):
                    content = f"Tool execution error: {result}"
                    logger.warning("Tool %s failed: %s", block.name, result)
                else:
                    content = result
                tool_results.append(
                    {
                        "type": "tool_result",
                        "tool_use_id": block.id,
                        "content": content,
                    }
                )

            messages.append({"role": "assistant", "content": response.content})
            messages.append({"role": "user", "content": tool_results})
            logger.info(
                "Agent iteration %d: called %s",
                iteration + 1,
                [b.name for b in tool_use_blocks],
            )
        else:
            # Unexpected stop reason — return whatever text is available
            logger.warning("Unexpected stop_reason: %s", response.stop_reason)
            return _extract_text(response) or f"Agent stopped unexpectedly: {response.stop_reason}"

    return (
        "The agent reached the maximum number of iterations without producing a final answer. "
        "Try a more specific question or call individual tools directly."
    )


async def _execute_tool(block: Any, registry: dict[str, ToolSpec]) -> str:
    """Look up a tool by name and execute it with the provided input kwargs."""
    spec = registry.get(block.name)
    if spec is None:
        return f"Unknown tool: '{block.name}'. Available tools: {list(registry.keys())}"
    try:
        return await spec.fn(**block.input)
    except Exception as exc:
        logger.warning("Tool '%s' raised: %s", block.name, exc)
        return f"Tool '{block.name}' error: {exc}"


def _extract_text(response: Any) -> str:
    """Extract concatenated text from a Claude response."""
    parts = [block.text for block in response.content if hasattr(block, "text")]
    return "\n".join(parts).strip()


def format_registry_docs(registry: dict[str, ToolSpec]) -> str:
    """Format the tool registry as a Markdown reference document.

    Intended for the ``tool://registry`` MCP resource so agents and clients
    can introspect available tools without invoking any of them.

    Tools are grouped by ``tool_category`` and sorted alphabetically within
    each group.  Each entry shows the tool name, description, ``use_when``
    guidance, and required inputs.

    Args:
        registry: Tool registry from ``build_tool_registry()``.

    Returns:
        Markdown string suitable for direct LLM or human consumption.
    """
    from collections import defaultdict

    _CATEGORY_ORDER = [
        "gene_annotation",
        "disease_evidence",
        "druggability",
        "structure",
        "pathways",
        "protein_engineering",
        "antibody",
        "synthesis",
    ]

    by_category: dict[str, list[ToolSpec]] = defaultdict(list)
    for spec in registry.values():
        by_category[spec.tool_category].append(spec)

    lines = [
        "# genesis-bio-mcp Tool Registry",
        "",
        f"{len(registry)} tools across {len(by_category)} categories.",
        "",
        "> **Semantic retrieval note:** The `use_when` field on each tool is written to be",
        "> embedding-searchable. At scale (50+ tools), embed all `use_when` strings once at",
        "> startup and retrieve the top-k most relevant tools per query via cosine similarity",
        "> before passing them to Claude. The `tool_category` field provides a coarse pre-filter.",
        "",
    ]

    ordered_categories = _CATEGORY_ORDER + sorted(
        c for c in by_category if c not in _CATEGORY_ORDER
    )

    for category in ordered_categories:
        specs = sorted(by_category.get(category, []), key=lambda s: s.name)
        if not specs:
            continue
        lines += [f"## {category.replace('_', ' ').title()} ({len(specs)} tools)", ""]
        for spec in specs:
            required = spec.input_schema.get("required", [])
            props = spec.input_schema.get("properties", {})
            param_strs = [f"`{p}` ({props[p].get('type', 'any')})" for p in required if p in props]
            params_summary = ", ".join(param_strs) if param_strs else "none"
            lines += [
                f"### `{spec.name}`",
                f"{spec.description}",
                f"- **Use when:** {spec.use_when}",
                f"- **Required inputs:** {params_summary}",
                "",
            ]

    return "\n".join(lines).rstrip()
