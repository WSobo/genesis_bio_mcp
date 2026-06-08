"""genesis_bio_mcp MCP server.

Exposes 44 tools for biomedical database queries:
  - resolve_gene                  UniProt + NCBI: gene symbol → canonical IDs
  - get_protein_info              UniProt Swiss-Prot protein annotation
  - get_protein_sequence          UniProt FASTA + biochem + liability scan
  - get_target_disease_association Open Targets: target–disease association score
  - get_cancer_dependency         DepMap: CRISPR essentiality across cancer lines
  - get_gwas_evidence             GWAS Catalog: genetic associations for a trait
  - get_compounds                 PubChem: active small molecules against a target
  - compute_molecular_properties  RDKit: physicochemical + drug-likeness from a SMILES (local)
  - predict_admet                 RDKit: ADMET liability panel (hERG/CYP/ESOL/alerts/SAscore) (local)
  - assess_selectivity            ChEMBL + corpus: measured + predicted off-target / selectivity profile
  - standardize_structure         RDKit: salt-strip/neutralize/canonical-tautomer + InChIKey (local)
  - search_similar_compounds      PubChem + RDKit: 2D-similarity / substructure structure search
  - get_chembl_compounds          ChEMBL: IC50/Ki/Kd + assay context (type, organism, confidence)
  - get_protein_structure         AlphaFold + RCSB PDB: structural data
  - get_structure_confidence      AlphaFold: per-residue pLDDT profile + disordered regions
  - get_structural_homologs       Foldseek: proteins with a similar 3D fold (structural search)
  - design_sequence_for_structure UMA-Inverse: redesign a backbone's sequence (inverse folding)
  - score_structure               UMA-Inverse: score a sequence vs a structure + mutation candidates
  - get_protein_interactome       STRING: binding partners and selectivity risks
  - get_biogrid_interactions      BioGRID: curated literature PPI network
  - get_antibody_structures       SAbDab: antibody/nanobody structures for an antigen
  - get_epitope_data              IEDB: known B-cell epitopes for an antigen
  - get_mhc_binding                IEDB NextGen: MHC-I/II binding prediction (T-cell epitopes)
  - get_cdr_developability        Antibody CDR liability scan (AbNum CDRs + biochem)
  - get_variant_constraints       gnomAD: gene-level LoF and missense constraint metrics
  - get_variant_effects           MyVariant + gnomAD + MaveDB + Ensembl VEP: pathogenicity + VEP
  - get_variant_consequences      Ensembl VEP: splice/UTR/regulatory + SIFT/PolyPhen
  - get_tissue_expression         GTEx: bulk-RNA median TPM across ~54 tissues
  - get_protein_atlas             HPA: tissue specificity, subcellular, pathology
  - get_domain_annotation         InterPro: domain boundaries, Pfam/SMART, GO terms
  - get_dms_scores                MaveDB: deep mutational scanning score sets
  - get_dms_variant_score         MaveDB: measured DMS score for one specific variant
  - get_drug_history              DGIdb + ClinicalTrials.gov + OpenFDA: drugs, trials, safety
  - get_pathway_context           Reactome: pathway membership and enrichment for a gene
  - get_pathway_members           Reactome: enumerate all genes in a named pathway
  - prioritize_target             Orchestration: full target assessment report
  - compare_targets               Compare 2–5 targets side by side for an indication
  - batch_resolve_genes           Resolve many gene names/aliases to canonical IDs in one call
  - batch_compute_molecular_properties RDKit: drug-likeness for many SMILES in one call (local)
  - corpus_describe                Corpus (Postgres+pgvector): manifest — coverage + provenance
  - corpus_search_targets_by_sequence Corpus: ESM-2 embedding kNN — targets similar to a query
  - corpus_find_similar_compounds  Corpus: Morgan/Tanimoto kNN — compounds similar to a SMILES
  - corpus_search_compounds        Corpus: HYBRID SQL filter + Tanimoto over bioactive compounds
  - run_biology_workflow          AI agent: dynamic tool selection for multi-step questions

All tools return Markdown strings for direct LLM consumption.
"""

from __future__ import annotations

import asyncio
import json as _json
import logging
import os
import time
from contextlib import asynccontextmanager
from datetime import UTC, datetime
from types import SimpleNamespace
from typing import Literal

import httpx
from mcp.server.fastmcp import FastMCP
from mcp.types import ToolAnnotations
from pydantic import BaseModel, ConfigDict, Field, model_validator

from genesis_bio_mcp import __version__
from genesis_bio_mcp.clients.alphafold import AlphaFoldClient
from genesis_bio_mcp.clients.biogrid import BioGRIDClient
from genesis_bio_mcp.clients.chembl import ChEMBLClient
from genesis_bio_mcp.clients.clinical_trials import ClinicalTrialsClient
from genesis_bio_mcp.clients.depmap import DepMapClient, load_depmap_cache
from genesis_bio_mcp.clients.dgidb import DGIdbClient
from genesis_bio_mcp.clients.ensembl import EnsemblClient
from genesis_bio_mcp.clients.foldseek import FoldseekClient
from genesis_bio_mcp.clients.gnomad import GnomADClient
from genesis_bio_mcp.clients.gtex import GTExClient
from genesis_bio_mcp.clients.gwas import GwasClient
from genesis_bio_mcp.clients.hpa import HPAClient
from genesis_bio_mcp.clients.iedb import IEDBClient
from genesis_bio_mcp.clients.iedb_tools import IEDBToolsClient
from genesis_bio_mcp.clients.interpro import InterProClient
from genesis_bio_mcp.clients.mavedb import MaveDBClient
from genesis_bio_mcp.clients.myvariant import MyVariantClient
from genesis_bio_mcp.clients.open_targets import OpenTargetsClient
from genesis_bio_mcp.clients.openfda import OpenFDAClient
from genesis_bio_mcp.clients.pubchem import PubChemClient
from genesis_bio_mcp.clients.reactome import ReactomeClient
from genesis_bio_mcp.clients.sabdab import SAbDabClient
from genesis_bio_mcp.clients.string_db import StringDbClient
from genesis_bio_mcp.clients.uma_inverse import UMAInverseClient
from genesis_bio_mcp.clients.uniprot import UniProtClient
from genesis_bio_mcp.clients.variant_effects import VariantEffectsClient
from genesis_bio_mcp.config.efo_resolver import EFOResolver
from genesis_bio_mcp.config.settings import settings
from genesis_bio_mcp.corpus import (
    create_corpus_pool,
    fetch_manifest,
    hybrid_search_compounds,
    resolve_target_accession,
    search_similar_targets,
)
from genesis_bio_mcp.corpus import (
    search_similar_compounds as _corpus_search_similar_compounds,
)
from genesis_bio_mcp.models import (
    BatchGeneResolution,
    BatchGeneResolutionItem,
    BatchMolecularProperties,
    ComparisonReport,
    CorpusCompoundHit,
    CorpusCompoundSearchHit,
    CorpusCompoundSearchResults,
    CorpusManifest,
    CorpusSimilarCompounds,
    CorpusTargetHit,
    CorpusTargetNeighbors,
    DMSVariantLookup,
    DrugHistory,
    ProteinSequence,
    Provenance,
    TargetComparisonRow,
    VariantEffects,
)
from genesis_bio_mcp.tools.admet import assess_admet as _assess_admet
from genesis_bio_mcp.tools.biochem import compute_features, scan_liabilities
from genesis_bio_mcp.tools.cdr_developability import assess_cdr_developability
from genesis_bio_mcp.tools.cheminformatics import (
    compute_molecular_properties as _compute_molecular_properties,
)
from genesis_bio_mcp.tools.cheminformatics import (
    morgan_fp_bits as _morgan_fp_bits,
)
from genesis_bio_mcp.tools.cheminformatics import (
    standardize_structure as _standardize_structure,
)
from genesis_bio_mcp.tools.gene_resolver import resolve_gene as _resolve_gene
from genesis_bio_mcp.tools.health import build_health_report, check_upstreams, health_to_markdown
from genesis_bio_mcp.tools.selectivity import assess_selectivity as _assess_selectivity
from genesis_bio_mcp.tools.target_prioritization import (
    attach_safety_signals as _attach_safety_signals,
)
from genesis_bio_mcp.tools.target_prioritization import (
    build_comparison_row as _build_comparison_row,
)
from genesis_bio_mcp.tools.target_prioritization import (
    prioritize_target as _prioritize_target,
)
from genesis_bio_mcp.tools.variant_parser import canonical_three_letter, parse_protein_change
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


def build_app_state(
    client: httpx.AsyncClient,
    gene_dep_cache: dict | None = None,
    *,
    corpus_pool: object | None = None,
) -> SimpleNamespace:
    """Construct the shared client graph that every tool reads from ``mcp.state``.

    Extracted from ``lifespan`` so the eval harness (and any out-of-server caller)
    can wire the exact same ~25 clients against a live ``httpx.AsyncClient`` without
    duplicating constructors. The async preludes — DepMap cache (``load_depmap_cache``)
    and corpus pool (``create_corpus_pool``) — are the caller's responsibility and
    passed in. ``gene_dep_cache`` defaults to empty (DepMap then uses its OT proxy).
    """
    state = SimpleNamespace()
    state.uniprot = UniProtClient(client)
    state.open_targets = OpenTargetsClient(client)
    state.depmap = DepMapClient(client, gene_dep_cache or {})
    state.gwas = GwasClient(client, efo_resolver=EFOResolver(client))
    state.gnomad = GnomADClient(client)
    state.interpro = InterProClient(client)
    state.pubchem = PubChemClient(client)
    state.chembl = ChEMBLClient(client)
    state.alphafold = AlphaFoldClient(client)
    state.foldseek = FoldseekClient(client, alphafold=state.alphafold)
    state.uma_inverse = UMAInverseClient(client)
    state.string_db = StringDbClient(client)
    state.biogrid = BioGRIDClient(client)
    state.sabdab = SAbDabClient(client)
    state.iedb = IEDBClient(client)
    state.iedb_tools = IEDBToolsClient(client)
    state.mavedb = MaveDBClient(client)
    state.myvariant = MyVariantClient(client)
    state.ensembl = EnsemblClient(client)
    state.gtex = GTExClient(client, ensembl=state.ensembl)
    state.hpa = HPAClient(client)
    state.variant_effects = VariantEffectsClient(
        gnomad=state.gnomad,
        myvariant=state.myvariant,
        mavedb=state.mavedb,
        ensembl=state.ensembl,
    )
    state.dgidb = DGIdbClient(client)
    state.clinical_trials = ClinicalTrialsClient(client)
    state.openfda = OpenFDAClient(client)
    state.reactome = ReactomeClient(client)
    state.corpus_pool = corpus_pool
    return state


@asynccontextmanager
async def lifespan(server: FastMCP):
    """Manage a shared httpx.AsyncClient and pre-load the DepMap gene cache."""
    if not os.environ.get("ANTHROPIC_API_KEY"):
        logger.warning(
            "ANTHROPIC_API_KEY is not set. run_biology_workflow will fail. "
            "Set it in claude_desktop_config.json under 'env' or export it before "
            "starting the MCP server."
        )

    async with httpx.AsyncClient(
        headers=_HEADERS,
        timeout=httpx.Timeout(settings.httpx_timeout, pool=settings.httpx_pool_timeout_secs),
        limits=httpx.Limits(
            max_connections=settings.httpx_max_connections,
            max_keepalive_connections=settings.httpx_max_keepalive,
        ),
        follow_redirects=True,
    ) as client:
        # Fetch DepMap gene_dep_summary once at startup for instant lookups
        gene_dep_cache = await load_depmap_cache(client)
        # Optional embedding-backed corpus store (v0.6.0). None when GENESIS_CORPUS_DSN is
        # unset or the DB is unreachable — corpus_* tools degrade gracefully, rest is unaffected.
        corpus_pool = await create_corpus_pool()

        server.state = build_app_state(client, gene_dep_cache, corpus_pool=corpus_pool)

        # Runtime metadata for the health://status resource.
        server.state.http_client = client
        server.state.started_at = time.time()
        server.state.depmap_gene_count = len(gene_dep_cache)
        try:
            yield
        finally:
            if server.state.corpus_pool is not None:
                await server.state.corpus_pool.close()


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
        resolution = await _resolve_gene(
            gene_name,
            uniprot_client=mcp.state.uniprot,
            ensembl_client=getattr(mcp.state, "ensembl", None),
        )
        symbol = resolution.hgnc_symbol or gene_name.strip().upper()
        return symbol, resolution.ncbi_gene_id
    except Exception:
        return gene_name.strip().upper(), None


# Provenance source labels keyed by output-model class name. Keyed by name (not the
# class object) so no extra imports are needed; an unmapped model falls back to the
# generic server label rather than failing — provenance degrades gracefully.
_SOURCE_BY_MODEL: dict[str, str] = {
    "GeneResolution": "UniProt + NCBI",
    "BatchGeneResolution": "UniProt + NCBI",
    "ProteinInfo": "UniProt",
    "ProteinSequence": "UniProt",
    "TargetDiseaseAssociation": "Open Targets",
    "CancerDependency": "DepMap",
    "GwasEvidence": "GWAS Catalog",
    "Compounds": "PubChem",
    "MolecularProperties": "RDKit (local)",
    "ADMETProfile": "RDKit (local)",
    "SelectivityProfile": "ChEMBL + corpus (RDKit similarity)",
    "BatchMolecularProperties": "RDKit (local)",
    "StandardizedStructure": "RDKit (local)",
    "SimilarCompounds": "PubChem + RDKit",
    "ChEMBLCompounds": "ChEMBL",
    "ProteinStructure": "AlphaFold + RCSB PDB",
    "StructureConfidence": "AlphaFold",
    "StructuralHomologs": "Foldseek",
    "SequenceDesign": "UMA-Inverse",
    "StructureScore": "UMA-Inverse",
    "ProteinInteractome": "STRING",
    "BioGRIDInteractome": "BioGRID",
    "EpitopeResults": "IEDB",
    "MHCBindingResults": "IEDB NextGen Tools",
    "AntibodyStructures": "SAbDab",
    "CDRDevelopabilityReport": "AbNum (Chothia) + biochem (local)",
    "GnomADConstraint": "gnomAD",
    "VariantEffects": "gnomAD + MyVariant.info + MaveDB + Ensembl VEP",
    "VEPConsequenceReport": "Ensembl VEP",
    "TissueExpressionProfile": "GTEx",
    "ProteinAtlasReport": "Human Protein Atlas",
    "DomainAnnotations": "InterPro",
    "DMSResults": "MaveDB",
    "DMSVariantLookup": "MaveDB",
    "DrugHistory": "DGIdb + ClinicalTrials.gov + OpenFDA",
    "PathwayContext": "Reactome",
    "TargetPrioritizationReport": "genesis-bio-mcp (multi-source synthesis)",
    "ComparisonReport": "genesis-bio-mcp (multi-source synthesis)",
    "CorpusManifest": "genesis-bio-mcp corpus (Postgres + pgvector)",
    "CorpusTargetNeighbors": "genesis-bio-mcp corpus (Postgres + pgvector)",
    "CorpusSimilarCompounds": "genesis-bio-mcp corpus (Postgres + pgvector)",
    "CorpusCompoundSearchResults": "genesis-bio-mcp corpus (Postgres + pgvector)",
}

# Attribute names a result model may expose for each provenance field, in priority
# order. The first present, non-empty value wins. Lets provenance pull the resolved
# query / upstream version / confidence straight off the model with no per-tool wiring.
_QUERY_ATTRS = (
    "gene_symbol",
    "hgnc_symbol",
    "query_gene",
    "query_smiles",
    "input_smiles",
    "antigen_query",
    "target",
)
_VERSION_ATTRS = ("dataset_version", "release", "data_version", "source_version", "model_version")
_CONFIDENCE_ATTRS = ("confidence", "mean_confidence")


def _first_attr(result: object, attrs: tuple[str, ...]) -> object | None:
    for attr in attrs:
        value = getattr(result, attr, None)
        if value not in (None, ""):
            return value
    return None


def _now_iso() -> str:
    """Current UTC time as an ISO-8601 'Z' timestamp (second precision)."""
    return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def _build_provenance(result: object) -> Provenance:
    """Derive the provenance block for a result from its type + common fields."""
    source = "genesis-bio-mcp"
    if result is not None:
        source = _SOURCE_BY_MODEL.get(type(result).__name__, "genesis-bio-mcp")
    query = _first_attr(result, _QUERY_ATTRS)
    version = _first_attr(result, _VERSION_ATTRS)
    confidence = _first_attr(result, _CONFIDENCE_ATTRS)
    return Provenance(
        source=source,
        source_version=str(version) if version is not None else None,
        retrieved_at=_now_iso(),
        query=str(query) if query is not None else None,
        confidence=float(confidence) if isinstance(confidence, (int, float)) else None,
    )


# Typed error taxonomy (agent contract) — machine-readable status an agent can
# branch on deterministically instead of string-matching prose. Surfaced in JSON
# error output only; Markdown error text is unchanged.
#   NotFound            — query was valid but no matching data exists upstream (default)
#   InvalidInput        — the input itself was malformed (bad SMILES, unparseable mutation)
#   RateLimited         — an upstream throttled the request (reserved; clients currently retry)
#   UpstreamUnavailable — an upstream service failed, timed out, or rejected the request
ErrorStatus = Literal["NotFound", "InvalidInput", "RateLimited", "UpstreamUnavailable"]


def _fmt(
    result: object,
    response_format: str,
    error_msg: str,
    *,
    error_status: ErrorStatus = "NotFound",
) -> str:
    """Format a Pydantic model as Markdown or JSON, or return an error representation.

    For ``response_format="json"`` the model is wrapped in a consistent agent contract
    envelope — ``{"provenance": {...}, "data": {...}}`` on success, or
    ``{"provenance": {...}, "error": {"status": ..., "message": ...}}`` when there is no
    result — so an agent always knows the source, timestamp, resolved query, and a
    machine-readable error status it can branch on. Markdown output is unchanged.

    Args:
        result: Pydantic model with .to_markdown() and .model_dump_json(), or None.
        response_format: "markdown" (default) or "json".
        error_msg: Human-readable error used when result is None.
        error_status: Typed status for the no-result case (see ``ErrorStatus``). Defaults
            to ``"NotFound"`` — the common "query valid, no data" case.
    """
    if result is None:
        if response_format == "json":
            return _json.dumps(
                {
                    "provenance": _build_provenance(None).model_dump(),
                    "error": {"status": error_status, "message": error_msg},
                },
                indent=2,
            )
        return f"**Error:** {error_msg}"
    if response_format == "json":
        envelope = {
            "provenance": _build_provenance(result).model_dump(),
            "data": _json.loads(result.model_dump_json()),  # type: ignore[attr-defined]
        }
        return _json.dumps(envelope, indent=2)
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


class GetProteinSequenceInput(_GeneInput):
    """Input for get_protein_sequence."""

    start: int | None = Field(
        default=None,
        description=(
            "Optional 1-indexed inclusive slice start. Use together with 'end' to "
            "request a specific domain or region (e.g. kinase domain 457-717)."
        ),
        ge=1,
    )
    end: int | None = Field(
        default=None,
        description="Optional 1-indexed inclusive slice end. Must be ≥ start.",
        ge=1,
    )

    @model_validator(mode="after")
    def _validate_region(self) -> GetProteinSequenceInput:
        if (self.start is None) != (self.end is None):
            raise ValueError("start and end must both be provided or both omitted")
        if self.start is not None and self.end is not None and self.end < self.start:
            raise ValueError("end must be >= start")
        return self


class GetCancerDependencyInput(_GeneInput):
    """Input for get_cancer_dependency."""


class GetCompoundsInput(_GeneInput):
    """Input for get_compounds."""


class ComputeMolecularPropertiesInput(BaseModel):
    """Input for compute_molecular_properties."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    smiles: str = Field(
        ...,
        description="SMILES string of the molecule, e.g. 'CC(=O)Oc1ccccc1C(=O)O' (aspirin).",
        min_length=1,
        max_length=2000,
    )
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD


class PredictADMETInput(BaseModel):
    """Input for predict_admet."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    smiles: str = Field(
        ...,
        description="SMILES string of the molecule, e.g. 'CC(=O)Oc1ccccc1C(=O)O' (aspirin).",
        min_length=1,
        max_length=2000,
    )
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD


class AssessSelectivityInput(BaseModel):
    """Input for assess_selectivity."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    smiles: str = Field(
        ...,
        description="SMILES of the compound to profile, e.g. a kinase inhibitor.",
        min_length=1,
        max_length=2000,
    )
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD


class StandardizeStructureInput(BaseModel):
    """Input for standardize_structure."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    smiles: str = Field(
        ...,
        description="SMILES string to standardize (salts/solvents stripped, charges neutralized).",
        min_length=1,
        max_length=2000,
    )
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD


class SearchSimilarCompoundsInput(BaseModel):
    """Input for search_similar_compounds."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    smiles: str = Field(
        ...,
        description="Query SMILES, e.g. 'CC(=O)Oc1ccccc1C(=O)O' (aspirin).",
        min_length=1,
        max_length=2000,
    )
    mode: Literal["similarity", "substructure"] = Field(
        "similarity",
        description="'similarity' = 2D Tanimoto ≥ threshold; 'substructure' = contains the query.",
    )
    threshold: float = Field(
        0.9,
        ge=0.4,
        le=1.0,
        description="Tanimoto similarity threshold (similarity mode only; 0.4–1.0).",
    )
    max_results: int = Field(20, ge=1, le=100, description="Maximum number of hits to return.")
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD


class BatchComputeMolecularPropertiesInput(BaseModel):
    """Input for batch_compute_molecular_properties."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    smiles_list: list[str] = Field(
        ...,
        description=(
            "List of SMILES strings to profile in one call (1–100). Unparseable entries "
            "are reported under 'failed' rather than failing the batch."
        ),
        min_length=1,
        max_length=100,
    )
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD


class BatchResolveGenesInput(BaseModel):
    """Input for batch_resolve_genes."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    gene_names: list[str] = Field(
        ...,
        description=(
            "List of gene symbols, aliases, or synonyms to resolve in one call (1–50). "
            "Examples: ['HER2', 'p53', 'COX2', 'ErbB1']."
        ),
        min_length=1,
        max_length=50,
    )
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD


class GetChEMBLCompoundsInput(_GeneInput):
    """Input for get_chembl_compounds."""


class GetProteinStructureInput(_GeneInput):
    """Input for get_protein_structure."""


class GetStructureConfidenceInput(_GeneInput):
    """Input for get_structure_confidence."""


class GetStructuralHomologsInput(_GeneInput):
    """Input for get_structural_homologs."""

    database: Literal["afdb-swissprot", "afdb-proteome", "pdb"] = Field(
        "afdb-swissprot",
        description=(
            "Foldseek target database: 'afdb-swissprot' (AlphaFold Swiss-Prot, curated, "
            "default) | 'afdb-proteome' (AlphaFold reference proteomes, broader) | 'pdb' "
            "(experimental PDB structures)."
        ),
    )
    max_hits: int | None = Field(
        None,
        ge=1,
        le=100,
        description="Maximum number of structural homologs to return. Omit for the server "
        "default (20).",
    )


_STRUCTURE_FIELD = Field(
    ...,
    description=(
        "The structure: full PDB text, OR a URL to a PDB file (e.g. the AlphaFold model URL "
        "reported by get_protein_structure, or an RCSB 'files.rcsb.org/download/{id}.pdb' URL). "
        "The live service caps inputs at 256 residues."
    ),
    min_length=1,
    max_length=2_000_000,
)


class DesignSequenceInput(BaseModel):
    """Input for design_sequence_for_structure."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    structure: str = _STRUCTURE_FIELD
    ligand: str | None = Field(
        None, description="Optional ligand context (3-letter code or SMILES).", max_length=2000
    )
    temperature: float = Field(
        0.1, ge=0.0, le=2.0, description="Sampling temperature; higher = more sequence diversity."
    )
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD


class ScoreStructureInput(BaseModel):
    """Input for score_structure."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    structure: str = _STRUCTURE_FIELD
    sequence: str | None = Field(
        None,
        description="Sequence to score against the structure (default: the structure's native "
        "sequence).",
        max_length=10000,
    )
    mode: Literal["autoregressive", "single-aa"] = Field(
        "autoregressive",
        description="'autoregressive' (fast; one pass) | 'single-aa' (per-residue; slower).",
    )
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD


class GetProteinInteractomeInput(_GeneInput):
    """Input for get_protein_interactome."""

    required_score: int | None = Field(
        None,
        ge=0,
        le=1000,
        description=(
            "STRING combined-score threshold 0–1000 for reported interactions. "
            "700 (default) = high confidence; lower it (e.g. 400 = medium) to surface "
            "more partners, raise it (e.g. 900 = highest) to keep only the strongest. "
            "Omit to use the server default."
        ),
    )
    limit: int | None = Field(
        None,
        gt=0,
        le=200,
        description=(
            "Maximum number of interaction partners to return. Omit to use the server default (20)."
        ),
    )


class GetBioGRIDInteractionsInput(_GeneInput):
    """Input for get_biogrid_interactions."""


class GetVariantConstraintsInput(_GeneInput):
    """Input for get_variant_constraints."""


class GetVariantEffectsInput(_GeneInput):
    """Input for get_variant_effects."""

    mutation: str = Field(
        ...,
        description=(
            "Protein change for the variant. Accepted forms: 'R175H', 'p.R175H', "
            "'Arg175His', 'p.Arg175His'. Case-insensitive; whitespace stripped."
        ),
        min_length=3,
        max_length=50,
    )


class GetDmsVariantScoreInput(_GeneInput):
    """Input for get_dms_variant_score."""

    mutation: str = Field(
        ...,
        description=(
            "Protein change to look up. Accepted forms: 'R175H', 'p.R175H', "
            "'Arg175His', 'p.Arg175His'. Case-insensitive; whitespace stripped."
        ),
        min_length=3,
        max_length=50,
    )


class GetVariantConsequencesInput(BaseModel):
    """Input for get_variant_consequences.

    Accepts any ONE of three query shapes (validated by the model_validator below):
      - ``gene_symbol`` + ``mutation`` (most common; resolves to canonical transcript HGVS.p)
      - ``hgvs_genomic`` (e.g. ``'7:g.140753336A>T'``)
      - ``chrom`` + ``pos`` + ``ref`` + ``alt``
    """

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")

    gene_symbol: str | None = Field(
        default=None,
        description="HGNC symbol paired with `mutation`. Required unless HGVS or coordinates are supplied.",
    )
    mutation: str | None = Field(
        default=None,
        description="Protein change paired with `gene_symbol`, e.g. 'V600E' or 'p.Val600Glu'.",
    )
    hgvs_genomic: str | None = Field(
        default=None,
        description="HGVS.g form, e.g. '7:g.140753336A>T'. Use when you already have coordinates.",
    )
    chrom: str | None = Field(
        default=None, description="Chromosome, e.g. '7' (paired with pos/ref/alt)."
    )
    pos: int | None = Field(
        default=None, description="1-indexed genomic position (paired with chrom/ref/alt)."
    )
    ref: str | None = Field(
        default=None, description="Reference allele (paired with chrom/pos/alt)."
    )
    alt: str | None = Field(
        default=None, description="Alternate allele (paired with chrom/pos/ref)."
    )
    include_all_transcripts: bool = Field(
        default=False,
        description="If True, return consequences across every transcript (can be 20+). Defaults to canonical only.",
    )
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD

    @model_validator(mode="after")
    def _exactly_one_query_shape(self) -> GetVariantConsequencesInput:
        shape_gene = bool(self.gene_symbol and self.mutation)
        shape_hgvs = bool(self.hgvs_genomic)
        shape_coord = all(v is not None for v in (self.chrom, self.pos, self.ref, self.alt))
        count = sum([shape_gene, shape_hgvs, shape_coord])
        if count != 1:
            raise ValueError(
                "Provide exactly one of: (gene_symbol + mutation), hgvs_genomic, "
                "or (chrom + pos + ref + alt)."
            )
        return self


class GetCDRDevelopabilityInput(BaseModel):
    """Input for get_cdr_developability.

    Provide at least one of:
      - ``vh`` and/or ``vl`` variable-domain sequences (auto-numbered to CDRs via AbNum, Chothia)
      - ``cdrs`` — explicit mapping of CDR field → sequence, e.g. ``{"vh_cdr3": "ARDYW"}``
    """

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")

    vh: str | None = Field(
        default=None,
        description="Heavy-chain variable domain (VH) sequence; auto-numbered to CDRs via AbNum (Chothia).",
    )
    vl: str | None = Field(
        default=None,
        description="Light-chain variable domain (VL) sequence; auto-numbered to CDRs via AbNum (Chothia).",
    )
    cdrs: dict[str, str] | None = Field(
        default=None,
        description="Explicit CDR sequences keyed by field: 'vh_cdr1','vh_cdr2','vh_cdr3','vl_cdr1','vl_cdr2','vl_cdr3'.",
    )
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD

    @model_validator(mode="after")
    def _at_least_one_input(self) -> GetCDRDevelopabilityInput:
        if not (self.vh or self.vl or self.cdrs):
            raise ValueError("Provide at least one of: vh, vl, or cdrs.")
        return self


class GetTissueExpressionInput(_GeneInput):
    """Input for get_tissue_expression."""


class GetProteinAtlasInput(_GeneInput):
    """Input for get_protein_atlas."""


class GetDomainAnnotationInput(_GeneInput):
    """Input for get_domain_annotation."""


class GetDMSScoresInput(_GeneInput):
    """Input for get_dms_scores."""


class GetMHCBindingInput(BaseModel):
    """Input for get_mhc_binding."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    sequence: str = Field(
        ...,
        description=(
            "Protein sequence, peptide, or multi-FASTA. Whole proteins are "
            "auto-windowed by the IEDB service into the requested peptide lengths. "
            "Examples: 'SLYNTVATL' (single 9mer), 'MAALSG...KLTQEHI' (protein region)."
        ),
        min_length=5,
        max_length=5000,
    )
    hla_alleles: list[str] | None = Field(
        default=None,
        description=(
            "HLA allele list, e.g. ['HLA-A*02:01', 'HLA-B*07:02']. If omitted, "
            "defaults to the 5-allele IEDB class-I reference panel (~85% global "
            "coverage) or class-II DRB1 panel (~60%) depending on mhc_class."
        ),
    )
    mhc_class: Literal["I", "II"] = Field(
        default="I",
        description="'I' for MHC Class I (CD8 / cytotoxic T cells) or 'II' for Class II (CD4 / helper).",
    )
    peptide_lengths: list[int] | None = Field(
        default=None,
        description=(
            "Peptide window lengths to test. Defaults: [9, 10] for class I, "
            "[15] for class II. Provide a single int or [min, max]."
        ),
    )
    method: str | None = Field(
        default=None,
        description=(
            "Predictor method. Default 'netmhcpan_el' (class I) or 'netmhciipan_el' (II). "
            "Other options: 'netmhcpan_ba' (binding affinity), 'mhcflurry', 'consensus'."
        ),
    )
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD


class GetEpitopeDataInput(BaseModel):
    """Input for get_epitope_data."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    antigen_query: str = Field(
        ...,
        description=(
            "Antigen protein name to search IEDB for known B-cell epitopes. "
            "Use the full protein name for best results: e.g. 'epidermal growth factor receptor', "
            "'tumor necrosis factor', 'programmed death-ligand 1'. "
            "Gene symbols (e.g. 'EGFR') also work but may return fewer results."
        ),
        min_length=1,
        max_length=200,
    )
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD


class GetAntibodyStructuresInput(BaseModel):
    """Input for get_antibody_structures."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    antigen_query: str = Field(
        ...,
        description=(
            "Antigen gene symbol or protein name to search for in SAbDab. "
            "Examples: 'EGFR', 'HER2', 'TNF', 'epidermal growth factor receptor', "
            "'programmed death-ligand 1', 'CD20'. Case-insensitive substring match."
        ),
        min_length=1,
        max_length=200,
    )
    max_results: int = Field(
        default=20,
        description="Maximum number of structures to return (default 20).",
        ge=1,
        le=100,
    )
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD


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

    Use this tool to understand a protein's biological function, Gene Ontology terms,
    enzyme class (EC), keywords, subcellular location, known disease-linked variants,
    and available 3D structures. Best used after resolve_gene to ensure a canonical
    symbol.

    Args:
        params (GetProteinInfoInput): gene_symbol, response_format.

    Returns:
        Markdown with function summary, sequence length, EC number, GO terms, keywords,
        pathways, disease associations, PDB IDs, known variants, and reviewed status.
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
async def get_protein_sequence(params: GetProteinSequenceInput) -> str:
    """Fetch the protein FASTA sequence and compute biochemistry + liability motifs.

    Call this tool FIRST for any sequence-level protein engineering question.
    Returns the amino-acid sequence (optionally sliced to a residue range),
    standard ExPASy ProtParam-style biochemistry (MW, pI, GRAVY, net charge,
    ε₂₈₀ reduced/oxidized), and a liability-motif scan covering:
    deamidation (NG/NS), isomerization (DG/DS), N-glycosylation sequons
    (N-X-S/T), oxidation-prone residues (M/W), and cysteines (flagged as
    free when UniProt annotates no disulfide bond at that position).

    Args:
        params (GetProteinSequenceInput): gene_symbol, optional start/end
            residue range, response_format.

    Returns:
        Markdown with biochemistry table, liability motifs grouped by type,
        annotated disulfide bonds, and a sequence preview.
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    protein = await mcp.state.uniprot.get_protein(symbol)
    if protein is None:
        return _fmt(
            None,
            params.response_format,
            f"No UniProt Swiss-Prot entry found for gene '{symbol}' in Homo sapiens.",
        )
    fetched = await mcp.state.uniprot.get_sequence(
        protein.uniprot_accession, start=params.start, end=params.end
    )
    if fetched is None:
        return _fmt(
            None,
            params.response_format,
            f"Could not fetch FASTA for '{symbol}' (UniProt {protein.uniprot_accession}).",
        )
    sequence, organism, description = fetched
    features = compute_features(sequence)
    # Translate full-sequence disulfide annotations into local positions if
    # a region was requested; otherwise use them directly.
    if params.start is not None and params.end is not None:
        offset = params.start - 1
        local_disulfide = {
            p - offset for p in protein.disulfide_bond_positions if params.start <= p <= params.end
        }
        region_disulfide_list = sorted(local_disulfide)
    else:
        local_disulfide = set(protein.disulfide_bond_positions)
        region_disulfide_list = list(protein.disulfide_bond_positions)
    liabilities = scan_liabilities(sequence, disulfide_annotated_positions=local_disulfide or None)
    result = ProteinSequence(
        uniprot_accession=protein.uniprot_accession,
        gene_symbol=protein.gene_symbol,
        organism=organism or protein.organism,
        description=description or protein.protein_name,
        sequence=sequence,
        region_start=params.start,
        region_end=params.end,
        features=features,
        liabilities=liabilities,
        disulfide_bond_positions=region_disulfide_list,
    )
    return _fmt(result, params.response_format, "")


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
        Markdown with overall_score and all per-datatype evidence scores
        (genetic_association, somatic_mutation, known_drug, literature,
        affected_pathway, rna_expression, animal_model), plus Open Targets
        tractability (druggability) buckets per modality (small molecule,
        antibody, PROTAC/degrader) — a target-level signal indicating how
        amenable the target is to each drug modality.
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

    Note: the DepMap summary endpoint reports per-gene dependency *counts*, not raw
    per-cell-line scores, so the reported mean dependency score is an approximated
    proxy (derived from the dependent-line fraction, or from OT somatic-mutation
    evidence on the fallback path) — NOT a measured CERES/Chronos value. The output
    labels this explicitly. The dependent-line fraction and pan-essential flag are
    the reliable signals.

    Args:
        params (GetCancerDependencyInput): gene_symbol, response_format.

    Returns:
        Markdown with fraction of dependent lines, pan-essential flag, top lineages,
        an explicitly-labeled approximated mean dependency score, and the data source
        (real DepMap counts or OT proxy).
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
        Markdown with total active compound count and top compounds by potency
        (IC50/EC50 in nM), each with name, molecular formula/weight, and the
        canonical SMILES structure string.
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    result = await mcp.state.pubchem.get_compounds(symbol)
    return _fmt(
        result, params.response_format, f"No PubChem bioactivity data found for gene '{symbol}'."
    )


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=False
    )
)
async def compute_molecular_properties(params: ComputeMolecularPropertiesInput) -> str:
    """Compute physicochemical and drug-likeness properties for a molecule from its SMILES.

    Locally computed with RDKit — deterministic, offline, no external API. Use this to
    profile a small molecule's developability without a database lookup: molecular weight,
    logP, TPSA, H-bond donors/acceptors, rotatable bonds, aromatic rings, fraction Csp3,
    QED, Lipinski Rule-of-Five and Veber checks, a PAINS assay-interference flag, and the
    Bemis-Murcko scaffold. Complements get_compounds / get_chembl_compounds (which retrieve
    known compounds) by analyzing any arbitrary structure.

    Args:
        params (ComputeMolecularPropertiesInput): smiles, response_format.

    Returns:
        Markdown with a property table, drug-likeness assessment (Lipinski/Veber/PAINS),
        and the Bemis-Murcko scaffold. Returns an error if the SMILES cannot be parsed.
    """
    result = _compute_molecular_properties(params.smiles)
    return _fmt(
        result,
        params.response_format,
        f"Invalid SMILES: '{params.smiles}' could not be parsed by RDKit.",
        error_status="InvalidInput",
    )


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=False
    )
)
async def predict_admet(params: PredictADMETInput) -> str:
    """Screen a molecule for ADMET / developability liabilities from its SMILES.

    Locally computed with RDKit — deterministic, offline, no external API. Returns a
    triage panel: Delaney **ESOL** aqueous solubility, GI-absorption / blood-brain-barrier
    heuristics, a **hERG** cardiotoxicity flag (basic amine + high lipophilicity), a
    **CYP3A4** metabolic-liability flag, Brenk + PAINS structural alerts, and Ertl **SAscore**
    synthetic accessibility. Use this to flag risks on a candidate or a virtual/enumerated
    structure before committing to it.

    IMPORTANT: every flag is a rule-based heuristic — a screening signal, NOT a validated
    prediction. Treat hERG/CYP/solubility calls as triage, not ground truth.

    Args:
        params (PredictADMETInput): smiles, response_format.

    Returns:
        Markdown panel of ADMET liability flags, structural alerts, and synthesizability.
        Returns an error if the SMILES cannot be parsed.
    """
    result = _assess_admet(params.smiles)
    return _fmt(
        result,
        params.response_format,
        f"Invalid SMILES: '{params.smiles}' could not be parsed by RDKit.",
        error_status="InvalidInput",
    )


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def assess_selectivity(params: AssessSelectivityInput) -> str:
    """Profile a compound's off-targets / kinome selectivity from its SMILES.

    Two honestly-separated layers:
    - **Measured** off-targets from ChEMBL (keyless) — every human target the compound has
      demonstrated pChEMBL bioactivity against, with a qualitative selectivity call
      (Selective / Moderately selective / Promiscuous).
    - **Predicted** off-target kinases from the optional kinome corpus — kinases hit by
      structurally similar known binders (chemical-similarity analogy), which also covers
      novel structures ChEMBL has never seen. Shown only when a corpus is configured.

    IMPORTANT: measured data demonstrates binding; **absence does not prove selectivity**
    (ChEMBL coverage is sparse). Predicted off-targets are a hypothesis, not a measurement.

    Args:
        params (AssessSelectivityInput): smiles, response_format.

    Returns:
        Markdown selectivity profile (measured + predicted off-targets). Error if the SMILES
        cannot be parsed.
    """
    result = await _assess_selectivity(
        params.smiles,
        chembl=mcp.state.chembl,
        corpus_pool=mcp.state.corpus_pool,
    )
    return _fmt(
        result,
        params.response_format,
        f"Invalid SMILES: '{params.smiles}' could not be parsed by RDKit.",
        error_status="InvalidInput",
    )


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=False
    )
)
async def standardize_structure(params: StandardizeStructureInput) -> str:
    """Standardize a molecule from its SMILES — salt/solvent strip, neutralize, canonical tautomer.

    Locally computed with RDKit (deterministic, offline). Produces a registration-ready
    canonical SMILES plus InChI / InChIKey for exact-structure matching and cross-source
    deduplication — the data-readiness primitive to run before comparing or joining compounds.

    Args:
        params (StandardizeStructureInput): smiles, response_format.

    Returns:
        Markdown with the standardized SMILES, InChIKey/InChI, molecular formula, and flags
        for whether a salt/solvent was removed or the structure otherwise changed. Returns an
        error if the SMILES cannot be parsed.
    """
    result = _standardize_structure(params.smiles)
    return _fmt(
        result,
        params.response_format,
        f"Invalid SMILES: '{params.smiles}' could not be parsed by RDKit.",
        error_status="InvalidInput",
    )


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def search_similar_compounds(params: SearchSimilarCompoundsInput) -> str:
    """Find compounds structurally related to a query molecule via PubChem.

    Two modes: 'similarity' returns compounds with 2D Tanimoto ≥ threshold (ranked by an
    RDKit Morgan Tanimoto computed against the query); 'substructure' returns compounds that
    contain the query as a substructure. Use this to expand chemical space around a hit,
    find analogs/me-too compounds, or scaffold-hop from a known structure.

    Args:
        params (SearchSimilarCompoundsInput): smiles, mode, threshold, max_results,
            response_format.

    Returns:
        Markdown table of hits (PubChem CID, name, formula, MW, Tanimoto). Returns a
        no-results message if nothing matches or the search could not complete.
    """
    result = await mcp.state.pubchem.search_similar(
        params.smiles,
        mode=params.mode,
        threshold=params.threshold,
        max_results=params.max_results,
    )
    if result is not None and result.total_found == 0:
        result = None
    return _fmt(
        result,
        params.response_format,
        f"No compounds found for '{params.smiles}' ({params.mode} search).",
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
        Markdown with ChEMBL target ID, best pChEMBL value, compound count, the
        most advanced clinical phase + mechanism(s) of action of drugs on this
        target, and a ranked table of top compounds with assay types and potency
        values.
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
async def get_structure_confidence(params: GetStructureConfidenceInput) -> str:
    """Retrieve the AlphaFold per-residue pLDDT confidence profile for a protein.

    Use this for protein engineering and construct design. It downloads the AlphaFold
    model and reports the fraction of residues in each confidence band plus a table of
    contiguous low-confidence (pLDDT < 70) stretches, which usually correspond to
    flexible or intrinsically disordered regions. Prefer well-ordered (pLDDT ≥ 70)
    regions when choosing mutation sites, truncation boundaries, or rigid scaffolds.
    Complements get_protein_structure, which reports only the mean pLDDT.

    Args:
        params (GetStructureConfidenceInput): gene_symbol, response_format.

    Returns:
        Markdown with mean pLDDT, per-band residue percentages, and a table of
        low-confidence regions.
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    protein = await mcp.state.uniprot.get_protein(symbol)
    accession = protein.uniprot_accession if protein else None
    result = await mcp.state.alphafold.get_confidence(symbol, uniprot_accession=accession)
    return _fmt(
        result,
        params.response_format,
        f"No AlphaFold per-residue confidence available for '{symbol}'. "
        "The protein may lack an AlphaFold model.",
    )


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_structural_homologs(params: GetStructuralHomologsInput) -> str:
    """Find proteins with a similar 3D fold via a Foldseek structural search.

    Submits the gene's AlphaFold model to the Foldseek web server and returns the
    most structurally-similar proteins. Unlike sequence search, this surfaces
    distant structural relatives whose sequence has diverged — useful for finding
    candidate scaffolds, inferring function for poorly-annotated targets, spotting
    structurally-similar off-target folds, and identifying fold families.

    Note: this runs an asynchronous remote search (submit → poll → fetch) and can
    take tens of seconds to a couple of minutes. The query's own AlphaFold model is
    excluded from the results.

    Args:
        params (GetStructuralHomologsInput): gene_symbol, database (afdb-swissprot |
            afdb-proteome | pdb), max_hits, response_format.

    Returns:
        Markdown table of structural homologs sorted by E-value: target id,
        description, E-value, bit score, match probability, and sequence identity.
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    protein = await mcp.state.uniprot.get_protein(symbol)
    accession = protein.uniprot_accession if protein else None
    result = await mcp.state.foldseek.search(
        symbol, accession, database=params.database, max_hits=params.max_hits
    )
    return _fmt(
        result,
        params.response_format,
        f"No structural homologs found for '{symbol}'. The protein may lack an "
        "AlphaFold model, or the Foldseek search did not complete.",
    )


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=False, openWorldHint=True
    )
)
async def design_sequence_for_structure(params: DesignSequenceInput) -> str:
    """Redesign a protein backbone's sequence via the UMA-Inverse inverse-folding model.

    Given a structure (PDB text or a URL to one — e.g. the AlphaFold model URL from
    get_protein_structure), returns amino-acid sequence(s) the model predicts will fold
    into that backbone, with per-residue and mean confidence. Use this for protein
    engineering / de novo sequence design against a known fold, or to propose a stabilized
    or solubilized variant of a target. Pairs with score_structure (score → redesign loop).

    Note: this calls a deployed model service; latency scales with structure length (the
    live endpoint caps inputs at 256 residues). Sampling is stochastic — raise temperature
    for more diversity.

    Args:
        params (DesignSequenceInput): structure (PDB text or URL), optional ligand,
            temperature, response_format.

    Returns:
        Markdown with the designed sequence(s), mean confidence, and low-confidence
        positions. Returns an error if the structure can't be read or the service rejects
        it (e.g. over the residue cap) / is unavailable.
    """
    pdb = await mcp.state.uma_inverse.resolve_pdb(params.structure)
    if not pdb:
        return _fmt(
            None,
            params.response_format,
            "Could not read the structure (empty input, or the PDB URL could not be fetched).",
            error_status="InvalidInput",
        )
    result = await mcp.state.uma_inverse.design(
        pdb, ligand=params.ligand, temperature=params.temperature
    )
    return _fmt(
        result,
        params.response_format,
        "UMA-Inverse design did not complete — the structure may exceed the service's "
        "residue cap (256 on the live endpoint), or the service is unavailable.",
        error_status="UpstreamUnavailable",
    )


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=False, openWorldHint=True
    )
)
async def score_structure(params: ScoreStructureInput) -> str:
    """Score a sequence against a structure with UMA-Inverse, flagging candidate mutations.

    Given a structure (PDB text or URL) and an optional sequence (defaults to the
    structure's native sequence), returns how well that sequence fits the backbone —
    overall perplexity (lower = better fit) and recovery, plus per-residue probabilities and
    the model's preferred residue at each position. Positions where the model prefers a
    different residue are surfaced as a **candidate-mutation table**. Use this to identify
    suboptimal residues to engineer, then feed the structure to design_sequence_for_structure.

    Note: this calls a deployed model service (live endpoint caps inputs at 256 residues).

    Args:
        params (ScoreStructureInput): structure (PDB text or URL), optional sequence, mode
            (autoregressive | single-aa), response_format.

    Returns:
        Markdown with perplexity/recovery and a candidate-mutation table. Returns an error if
        the structure can't be read or the service rejects it / is unavailable.
    """
    pdb = await mcp.state.uma_inverse.resolve_pdb(params.structure)
    if not pdb:
        return _fmt(
            None,
            params.response_format,
            "Could not read the structure (empty input, or the PDB URL could not be fetched).",
            error_status="InvalidInput",
        )
    result = await mcp.state.uma_inverse.score(pdb, sequence=params.sequence, mode=params.mode)
    return _fmt(
        result,
        params.response_format,
        "UMA-Inverse scoring did not complete — the structure may exceed the service's "
        "residue cap (256 on the live endpoint), or the service is unavailable.",
        error_status="UpstreamUnavailable",
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
        params (GetProteinInteractomeInput): gene_symbol, response_format, and optional
            required_score (STRING confidence threshold 0–1000) and limit (max partners).

    Returns:
        Markdown with top interaction partners sorted by STRING confidence score,
        each partner's dominant evidence channel and its sub-score (so experimental /
        database support can be told apart from textmining co-mentions at the same
        combined score), the contributing evidence channels, and total partner count.
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    result = await mcp.state.string_db.get_interactome(
        symbol, required_score=params.required_score, limit=params.limit
    )
    if result is not None and result.total_partners == 0:
        result = None
    return _fmt(result, params.response_format, f"No STRING interaction data found for '{symbol}'.")


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_biogrid_interactions(params: GetBioGRIDInteractionsInput) -> str:
    """Retrieve curated protein–protein interactions from BioGRID.

    Use this tool to get literature-curated PPI data from BioGRID, which provides
    manually curated interaction records with experimental method metadata (two-hybrid,
    co-IP, proximity ligation, etc.) and PubMed citations.  Complements STRING (which
    scores confidence across multiple evidence types) with individual experiment-level
    records.

    Requires the BIOGRID_ACCESS_KEY environment variable.  Register for a free key at
    https://webservice.thebiogrid.org/

    Args:
        params (GetBioGRIDInteractionsInput): gene_symbol, response_format.

    Returns:
        Markdown table of top interaction partners ranked by evidence count, with
        experimental method and interaction type annotations.  Returns an error
        message if the API key is not set.
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    if not __import__("os").environ.get("BIOGRID_ACCESS_KEY"):
        return (
            "**BioGRID data unavailable:** `BIOGRID_ACCESS_KEY` is not set.\n\n"
            "Register for a free key at https://webservice.thebiogrid.org/ and add it to the "
            "server environment as `BIOGRID_ACCESS_KEY=<your-key>`."
        )
    result = await mcp.state.biogrid.get_interactions(symbol)
    return _fmt(
        result,
        params.response_format,
        f"No BioGRID interaction data found for '{symbol}'.",
    )


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_antibody_structures(params: GetAntibodyStructuresInput) -> str:
    """Retrieve antibody and nanobody (VHH) structures from SAbDab for a given antigen.

    Use this tool when designing or evaluating antibody or nanobody therapeutics, or when
    you need PDB-curated structural templates of antibodies bound to a target antigen.
    SAbDab provides CDR-annotated structures for both conventional Fabs/IgGs and VHH
    single-domain nanobodies.

    Searches antigen_name, compound, and antigen description fields (case-insensitive
    substring match).  No API key required; results are cached locally for 7 days.

    Args:
        params (GetAntibodyStructuresInput): antigen_query, max_results, response_format.

    Returns:
        Markdown table of matching antibody/nanobody structures sorted by resolution
        (best first), including PDB ID, type (Fab/VHH/scFv), resolution, experimental
        method, antibody species, germline subclass, engineered flag, and measured affinity
        where available.
    """
    result = await mcp.state.sabdab.get_antibody_structures(
        params.antigen_query, max_results=params.max_results
    )
    return _fmt(
        result,
        params.response_format,
        f"SAbDab data temporarily unavailable for '{params.antigen_query}'. "
        "The SAbDab summary TSV could not be downloaded. Try again later.",
        error_status="UpstreamUnavailable",
    )


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_epitope_data(params: GetEpitopeDataInput) -> str:
    """Retrieve known B-cell epitope records from IEDB for a given antigen.

    Use this tool to understand the immunological landscape for an antigen target
    before designing antibody or nanobody therapeutics.  IEDB provides the most
    comprehensive curated database of published B-cell epitopes with experimental
    evidence (binding assays, crystal structures).

    Critical for antibody design: known epitopes identify which regions of the
    antigen surface have been targeted, and PDB-annotated epitopes confirm the
    3D binding mode.  A target with many known epitopes (especially structural ones)
    is well-characterized for Ab engineering; a sparse epitope landscape suggests
    fewer published precedents.

    Use the full protein name for best results (e.g. 'epidermal growth factor
    receptor' not 'EGFR'). Combine with get_antibody_structures for structural
    context and get_protein_structure for antigen conformation.

    Args:
        params (GetEpitopeDataInput): antigen_query, response_format.

    Returns:
        Markdown with total positive assays, unique epitope count, structural
        evidence count, and a table of epitope sequences/residues with isotype,
        position, PDB ID, and PubMed citation.
    """
    result = await mcp.state.iedb.get_epitopes(params.antigen_query)
    return _fmt(
        result,
        params.response_format,
        f"No B-cell epitope data found for '{params.antigen_query}' in IEDB.",
    )


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_mhc_binding(params: GetMHCBindingInput) -> str:
    """Predict MHC-I or MHC-II binding for peptides via IEDB NextGen Tools.

    Submits a NetMHCpan (or user-chosen) prediction job to IEDB's async
    NextGen Tools pipeline, polls for completion, and returns per-peptide,
    per-allele predicted percentile ranks and binder classification.

    Use this tool for T-cell epitope / immunogenicity assessment:
      - Therapeutic antibody CDR3 → screen against human HLA panel to flag
        potential immunogenicity risk.
      - Vaccine antigen → identify candidate T-cell epitopes.
      - Engineered construct → predict presentation before advancing to
        assays.

    Percentile-rank interpretation (IEDB convention):
      * < 0.5%  strong binder, high confidence
      * < 2%    weak binder, likely to bind
      * ≥ 2%    non-binder

    This call can take up to 60 seconds due to the IEDB async job queue.
    Large requests (many peptides × many alleles) may time out with a note
    in the output — narrow the input and retry.

    Args:
        params (GetMHCBindingInput): sequence, hla_alleles, mhc_class,
            peptide_lengths, method, response_format.

    Returns:
        Markdown with binder counts, top-10 hits table sorted by percentile,
        and per-allele strong-binder distribution.
    """
    try:
        result = await mcp.state.iedb_tools.predict_mhc_binding(
            sequence=params.sequence,
            alleles=params.hla_alleles,
            mhc_class=params.mhc_class,
            peptide_lengths=params.peptide_lengths,
            method=params.method,
        )
    except ValueError as exc:
        return _fmt(None, params.response_format, str(exc), error_status="InvalidInput")
    return _fmt(
        result,
        params.response_format,
        "IEDB NextGen Tools prediction failed — the service may be "
        "temporarily unavailable. Retry shortly.",
        error_status="UpstreamUnavailable",
    )


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_variant_constraints(params: GetVariantConstraintsInput) -> str:
    """Retrieve gene-level evolutionary constraint metrics from gnomAD v4.

    Use this tool to understand how much mutation a gene tolerates in the human
    population before designing variants or selecting engineering targets.  High
    constraint (pLI > 0.9, LOEUF < 0.35) means the gene does not tolerate
    loss-of-function — avoid broad mutagenesis and prefer conservative substitutions
    in well-characterized tolerant regions.  Unconstrained genes (low pLI, oe_mis ≈ 1)
    support bolder engineering strategies.

    This is a critical pre-filter for any protein engineering campaign: run it before
    DMS lookup, variant design, or stability engineering.

    Args:
        params (GetVariantConstraintsInput): gene_symbol, response_format.

    Returns:
        Markdown with pLI, LOEUF, oe_lof, oe_mis, Z-scores, and an
        interpretation of engineering implications.
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    result = await mcp.state.gnomad.get_constraint(symbol)
    return _fmt(result, params.response_format, f"No gnomAD constraint data found for '{symbol}'.")


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_variant_effects(params: GetVariantEffectsInput) -> str:
    """Fetch mutation-level pathogenicity and fitness for a specific missense variant.

    Fans out to three sources in parallel:
      - **gnomAD v4** resolves the protein change to a canonical variant_id and
        returns population allele frequency.
      - **MyVariant.info** returns ClinVar submissions (clinical significance,
        review stars, conditions), gnomAD exome + genome allele frequency, the
        AlphaMissense (am_pathogenicity + class), REVEL, CADD, SIFT, PolyPhen-2,
        ESM1b, and EVE pathogenicity predictors, plus GERP++/phyloP evolutionary
        conservation.
      - **MaveDB** probes the top 3 DMS score sets for the gene for a
        per-variant fitness score.

    Use this tool when a specific residue change matters — e.g. assessing
    whether an engineered mutation is likely to break function, checking the
    clinical status of a candidate driver, or cross-referencing DMS fitness
    against population evidence.

    Args:
        params (GetVariantEffectsInput): gene_symbol, mutation (e.g. 'R175H'
            or 'p.Arg175His'), response_format.

    Returns:
        Markdown synthesis with ClinVar verdict, AlphaMissense class + score,
        REVEL/CADD/SIFT/PolyPhen, gnomAD allele frequency, and MaveDB DMS
        scores for this exact mutation.
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    try:
        result: VariantEffects = await mcp.state.variant_effects.get_effects(
            symbol, params.mutation
        )
    except ValueError as exc:
        return _fmt(None, params.response_format, str(exc), error_status="InvalidInput")
    return _fmt(result, params.response_format, "")


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_variant_consequences(params: GetVariantConsequencesInput) -> str:
    """Predict functional consequences of a variant via Ensembl VEP.

    Distinct from ``get_variant_effects`` (pathogenicity / clinical) — this
    tool answers *what does the variant hit*: splice donor/acceptor, 5'/3'
    UTR, regulatory region, coding change, and which transcripts are
    affected. Returns VEP's SIFT and PolyPhen predictions, which are
    complementary to the dbNSFP-sourced SIFT/PolyPhen surfaced by MyVariant.
    Reports regulatory feature overlap (promoter, enhancer, CTCF site) when
    the variant falls inside one.

    Use this when the *location* of a variant matters — splice-impacting
    intronic variants, UTR variants with regulatory implications, non-coding
    variants pulled from GWAS fine-mapping, or when dissecting transcript-
    isoform-specific effects.

    Args:
        params (GetVariantConsequencesInput): provide exactly one of:
            gene_symbol + mutation, hgvs_genomic, or chrom+pos+ref+alt.
            ``include_all_transcripts=True`` returns every transcript
            (default canonical only).

    Returns:
        Markdown with most-severe consequence, per-transcript effects,
        SIFT/PolyPhen predictions, and regulatory feature overlap.
    """
    ensembl = mcp.state.ensembl
    if params.gene_symbol and params.mutation:
        symbol, _ = await _resolve_symbol(params.gene_symbol)
        result = await ensembl.get_vep_consequences(
            symbol,
            params.mutation,
            include_all_transcripts=params.include_all_transcripts,
        )
        err = f"No VEP consequences found for {symbol} {params.mutation}."
    elif params.hgvs_genomic:
        result = await ensembl.get_vep_by_hgvs(
            params.hgvs_genomic,
            include_all_transcripts=params.include_all_transcripts,
        )
        err = f"No VEP consequences found for '{params.hgvs_genomic}'."
    else:
        # Ensembl region form: "{chrom}:{pos}-{pos}:{strand}". Strand 1 is
        # a safe default since VEP normalizes to forward-strand alleles.
        region = f"{params.chrom}:{params.pos}-{params.pos}:1"
        result = await ensembl.get_vep_by_region(
            region,
            params.alt or "",
            include_all_transcripts=params.include_all_transcripts,
        )
        err = (
            f"No VEP consequences found for {params.chrom}:{params.pos} {params.ref}>{params.alt}."
        )
    return _fmt(result, params.response_format, err)


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_tissue_expression(params: GetTissueExpressionInput) -> str:
    """Retrieve GTEx bulk-RNA median tissue expression (TPM) for a gene.

    Returns per-tissue median TPM across ~54 GTEx tissues (brain subregions,
    heart, liver, kidney, etc.), sorted by expression level. Use this to
    answer:
      - Is the gene expressed in the tissue where the disease manifests?
      - Is expression restricted enough to give a therapeutic window?
      - Which tissues might see on-target toxicity?

    GTEx uses GENCODE IDs; we resolve the HGNC symbol → Ensembl ID
    automatically before querying. Genes without a GENCODE mapping (some
    non-coding or retired annotations) return an empty profile.

    Args:
        params (GetTissueExpressionInput): gene_symbol, response_format.

    Returns:
        Markdown table with the top tissues by median TPM, led by an explicit
        top-tissue callout and a GTEx-native tissue-specificity flag (the top
        tissue's fold-enrichment over the median tissue: 'Tissue-restricted' /
        'Tissue-enhanced' / 'Broadly expressed') — restricted expression generally
        implies a wider therapeutic window.
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    result = await mcp.state.gtex.get_expression(symbol)
    return _fmt(result, params.response_format, f"No GTEx expression data found for '{symbol}'.")


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_protein_atlas(params: GetProteinAtlasInput) -> str:
    """Retrieve Human Protein Atlas expression, subcellular, and pathology data.

    Returns HPA's RNA tissue-specificity category (Tissue enriched /
    Group enriched / Tissue enhanced / Low tissue specificity / Not detected),
    the numerical specificity score, subcellular localization (IHC-based),
    the protein-class annotation with a derived membrane/secreted call
    (target accessibility — whether the protein is reachable by antibodies or
    other biologics), and prognostic cancer-outcome data per indication.

    Complementary to get_tissue_expression (GTEx): HPA adds protein-level
    localization, druggability/accessibility class, and pathology context
    that bulk RNA cannot provide.

    Args:
        params (GetProteinAtlasInput): gene_symbol, response_format.

    Returns:
        Markdown with specificity category, subcellular locations, protein
        class + membrane/secreted accessibility, enhanced tissues, and
        pathology cancer rows.
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    result = await mcp.state.hpa.get_report(symbol)
    return _fmt(
        result, params.response_format, f"No Human Protein Atlas data found for '{symbol}'."
    )


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_domain_annotation(params: GetDomainAnnotationInput) -> str:
    """Retrieve InterPro domain annotations for a protein.

    Use this tool to understand which domains a protein contains before planning
    any engineering campaign — mutating within a conserved domain core often
    disrupts function.  Returns all InterPro-integrated domain entries (Pfam,
    SMART, CDD, etc.) with residue positions, member database cross-references,
    and GO term annotations.

    This is complementary to get_protein_structure (3D context) and
    get_variant_constraints (mutation tolerance): domain boundaries define *where*
    engineering is possible; constraint data defines *how much* tolerance exists.

    Args:
        params (GetDomainAnnotationInput): gene_symbol, response_format.

    Returns:
        Markdown table of domains sorted by position, with InterPro accession,
        type, residue range, Pfam/SMART cross-references, and GO terms.
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    # Need the UniProt accession to query InterPro
    protein = await mcp.state.uniprot.get_protein(symbol)
    accession = protein.uniprot_accession if protein else symbol
    result = await mcp.state.interpro.get_domains(symbol, accession)
    return _fmt(result, params.response_format, f"No InterPro domain data found for '{symbol}'.")


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_dms_scores(params: GetDMSScoresInput) -> str:
    """Retrieve deep mutational scanning (DMS) score sets from MaveDB for a gene.

    Use this tool when you need residue-level fitness or function scores for a
    protein.  DMS datasets provide the highest-resolution mutation-tolerance signal
    available — every possible single amino acid change scored in a defined assay
    (binding, stability, growth, etc.).  When a DMS dataset exists for a target it
    should be consulted before proposing specific mutations in an engineering campaign.

    Not all genes have DMS data.  Returns an empty result (not an error) when no
    score sets are found — the absence is itself informative.

    Args:
        params (GetDMSScoresInput): gene_symbol, response_format.

    Returns:
        Markdown table of available DMS score sets with URNs, variant counts,
        UniProt accessions, publication dates, and PubMed IDs.  Returns a message
        indicating no data when MaveDB has no score sets for this gene.
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    result = await mcp.state.mavedb.get_dms_scores(symbol)
    return _fmt(
        result, params.response_format, f"No DMS score sets found in MaveDB for '{symbol}'."
    )


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_dms_variant_score(params: GetDmsVariantScoreInput) -> str:
    """Look up the measured deep-mutational-scanning (DMS) score for one specific variant from MaveDB.

    Use this when you have a specific mutation and want the experimentally measured
    functional/fitness score — the highest-resolution residue-level signal available.
    This is more targeted than get_dms_scores (which lists the score-set catalog for a
    gene) and lighter than get_variant_effects (which bundles ClinVar/AlphaMissense/
    gnomAD/VEP around the DMS score). Probes the gene's largest MaveDB score sets and
    returns the score from each one that measured this variant.

    Args:
        params (GetDmsVariantScoreInput): gene_symbol, mutation, response_format.

    Returns:
        Markdown with the canonical HGVS p. form and a table of measured DMS scores
        (one row per score set), or a message if no DMS data covers the variant.
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    try:
        orig, pos, new = parse_protein_change(params.mutation)
    except ValueError as exc:
        return _fmt(
            None,
            params.response_format,
            f"Could not parse mutation '{params.mutation}': {exc}",
            error_status="InvalidInput",
        )
    hgvs_p = canonical_three_letter(orig, pos, new)
    dms = await mcp.state.mavedb.get_dms_scores(symbol)
    scores = await mcp.state.mavedb.get_variant_scores_for_gene(symbol, hgvs_p)
    result = DMSVariantLookup(
        gene_symbol=symbol,
        mutation=params.mutation,
        canonical_hgvs_protein=hgvs_p,
        score_sets_available=dms.total_score_sets if dms else 0,
        scores=scores,
    )
    return _fmt(result, params.response_format, "")


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_cdr_developability(params: GetCDRDevelopabilityInput) -> str:
    """Assess antibody CDR developability liabilities.

    Use this when engineering or triaging an antibody/nanobody. Provide VH and/or VL
    variable-domain sequences (auto-numbered to the six CDRs via AbNum, Chothia scheme)
    and/or explicit CDR sequences. For each CDR it reports net charge, hydrophobicity
    (GRAVY), and sequence-motif liabilities (deamidation NG/NS, isomerization DG/DS,
    N-glycosylation sequons, Met/Trp oxidation, cysteines), then flags the high-priority
    developability risks (CDR motif chemistry, hydrophobic/long CDR-H3). Complements
    get_antibody_structures (structural templates) and get_mhc_binding (immunogenicity).

    Args:
        params (GetCDRDevelopabilityInput): vh and/or vl sequences, or explicit cdrs;
            response_format.

    Returns:
        Markdown with a per-CDR table and a list of flagged developability liabilities.
    """
    cdr_map: dict[str, str] = {}
    sources: list[str] = []
    if params.vh or params.vl:
        numbered = await mcp.state.sabdab.number_chains(params.vh, params.vl)
        auto = {k: v for k, v in numbered.items() if v}
        if auto:
            cdr_map.update(auto)
            sources.append("AbNum (Chothia)")
    if params.cdrs:
        explicit = {k: v for k, v in params.cdrs.items() if v}
        if explicit:
            cdr_map.update(explicit)
            sources.append("user-provided")
    if not cdr_map:
        return (
            "No CDR sequences could be determined. If you supplied VH/VL, AbNum numbering "
            "may have failed or timed out — provide explicit CDRs via `cdrs` instead."
        )
    report = assess_cdr_developability(cdr_map, numbering_source=" + ".join(sources))
    return _fmt(report, params.response_format, "")


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_drug_history(params: GetDrugHistoryInput) -> str:
    """Retrieve the drug development history and post-market safety profile.

    Use this tool to understand the clinical precedent for targeting a gene:
    what drugs already exist (approved or investigational), what indications have
    been pursued, how many clinical trials have targeted this gene, and — for
    approved drugs — the FAERS adverse-event profile, FDA boxed warnings, and
    recent recalls. Essential for first-in-class vs. best-in-class strategy and
    for anticipating class-level safety liabilities.

    Args:
        params (GetDrugHistoryInput): gene_symbol, response_format.

    Returns:
        Markdown with known drugs (type, phase, approval status, DGIdb interaction
        score), trial counts by phase, a table of recent clinical trials from
        ClinicalTrials.gov (with lead sponsor, interventions, and enrollment), and a
        safety-signals section listing boxed warnings and top adverse events for
        the highest-phase approved drugs (up to 5).
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    drugs, ct_result = await asyncio.gather(
        mcp.state.dgidb.get_drug_interactions(symbol),
        mcp.state.clinical_trials.get_trials(symbol),
    )
    ct_trials, ct_counts = ct_result
    drugs = await _attach_safety_signals(drugs, openfda=mcp.state.openfda)
    result = DrugHistory(
        gene_symbol=symbol,
        known_drugs=drugs,
        approved_drug_count=sum(1 for d in drugs if d.approved),
        trial_counts_by_phase=ct_counts,
        recent_trials=ct_trials[:10],
    )
    # Empty-but-valid result gets a friendly Markdown note; JSON still returns the
    # (empty) model wrapped in the provenance envelope so the contract stays uniform.
    if not drugs and not ct_trials and params.response_format != "json":
        return (
            f"**No drug history found for '{symbol}'.** This may be a first-in-class opportunity."
        )
    return _fmt(result, params.response_format, "")


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
        Markdown with top enriched Reactome pathways, enrichment p-values and
        Benjamini-Hochberg FDR (multiple-testing-corrected; FDR < 0.05 flagged
        as significant), pathway categories, and gene counts per pathway.
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
    PubChem, ChEMBL) and returns a structured Markdown report with a per-axis
    evidence profile — each axis (disease association, clinical validation, genetics,
    cancer dependency, chemical tractability, annotation) rated independently as
    strong/moderate/weak/none/n-a. There is deliberately no single composite score,
    so the cross-axis trade-off stays visible for go/no-go reasoning on a
    target–indication pair.

    Set extended=True for a deeper report that also includes AlphaFold structure data,
    STRING interactome (selectivity risks), drug history (DGIdb + ClinicalTrials.gov),
    Reactome pathway context, and an HPA tissue-specificity axis. Extended mode adds
    ~10–20s to query time.

    Args:
        params (PrioritizeTargetInput): gene_symbol, indication, extended, response_format.

    Returns:
        Markdown report with a per-axis evidence profile, a plain-language evidence
        summary, data-source coverage, and data gaps. With extended=True, also includes
        structure, interactors, drugs, pathways, and tissue specificity.
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
        openfda=mcp.state.openfda if params.extended else None,
        reactome=mcp.state.reactome if params.extended else None,
        hpa=mcp.state.hpa if params.extended else None,
    )
    return _fmt(result, params.response_format, "")


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=False, openWorldHint=True
    )
)
async def compare_targets(params: CompareTargetsInput) -> str:
    """Compare 2–5 drug targets side by side for a given therapeutic indication.

    Runs a full prioritize_target assessment for each gene in parallel and returns a
    Markdown evidence matrix (genes × axes, each cell rated strong/moderate/weak), a
    raw-signals table, and per-gene evidence summaries. Rows are ordered by strength of
    evidence (count of strong, then moderate axes) — a transparent ordering for
    convenience, not a validated composite score. Use this when choosing between
    multiple candidate targets for an indication.

    Args:
        params (CompareTargetsInput): gene_symbols (2–5), indication, response_format.

    Returns:
        Markdown evidence matrix + raw signals + per-gene summaries, ordered by
        strength of evidence.
    """
    gene_symbols = list(params.gene_symbols)
    dropped: list[str] = []
    if len(gene_symbols) > 5:
        dropped = gene_symbols[5:]
        gene_symbols = gene_symbols[:5]
        logger.warning("compare_targets: capped to 5 genes, dropped: %s", dropped)

    # Bound gene-level concurrency so 5 genes don't oversubscribe the small
    # per-source semaphores (ChEMBL=2, OT=3) and starve each other, and cap each
    # gene's wall-clock so one slow gene becomes a timed-out row instead of
    # wedging the whole comparison. A gene exceeding its budget raises
    # TimeoutError, which _build_comparison_row renders as a failed row.
    gene_sem = asyncio.Semaphore(settings.compare_gene_concurrency)

    async def _assess(sym: str):
        async with gene_sem:
            return await asyncio.wait_for(
                _prioritize_target(
                    gene_symbol=sym,
                    indication=params.indication,
                    uniprot=mcp.state.uniprot,
                    open_targets=mcp.state.open_targets,
                    depmap=mcp.state.depmap,
                    gwas=mcp.state.gwas,
                    pubchem=mcp.state.pubchem,
                    chembl=mcp.state.chembl,
                ),
                timeout=settings.prioritize_gene_budget_secs,
            )

    reports = await asyncio.gather(
        *[_assess(sym) for sym in gene_symbols],
        return_exceptions=True,
    )

    rows: list[TargetComparisonRow] = []
    for sym, report in zip(gene_symbols, reports):
        if isinstance(report, Exception):
            logger.warning("compare_targets: failed for %s — %s", sym, report)
        rows.append(_build_comparison_row(sym, report))

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
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def batch_resolve_genes(params: BatchResolveGenesInput) -> str:
    """Resolve many gene names, aliases, or synonyms to canonical identifiers in one call.

    Agent fan-out helper: instead of calling resolve_gene N times, pass a whole list and
    get back one table mapping each query to its HGNC symbol, NCBI Gene ID, and UniProt
    accession. Resolutions run concurrently; an unresolved entry is flagged in its row
    rather than failing the batch. Use this to normalize a gene list before querying other
    tools.

    Args:
        params (BatchResolveGenesInput): gene_names (1–50), response_format.

    Returns:
        Markdown table with one row per query (symbol, NCBI Gene, UniProt, resolved status).
    """

    async def _one(name: str) -> BatchGeneResolutionItem:
        try:
            res = await _resolve_gene(name, uniprot_client=mcp.state.uniprot)
        except Exception:
            res = None
        resolved = bool(
            res
            and (res.source != "input" or res.ncbi_gene_id or res.uniprot_accession or res.hgnc_id)
        )
        return BatchGeneResolutionItem(query=name, resolved=resolved, resolution=res)

    items = await asyncio.gather(*[_one(n) for n in params.gene_names])
    batch = BatchGeneResolution(items=list(items), total_requested=len(params.gene_names))
    return _fmt(batch, params.response_format, "")


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=False
    )
)
async def batch_compute_molecular_properties(
    params: BatchComputeMolecularPropertiesInput,
) -> str:
    """Compute drug-likeness properties for many molecules from their SMILES in one call.

    Agent fan-out helper for the local RDKit compute_molecular_properties tool: pass a list
    of SMILES and get one compact summary table (formula, MW, logP, TPSA, QED, Lipinski Ro5,
    PAINS) with one row per molecule. Deterministic, offline, no external API. Unparseable
    SMILES are listed under 'failed' instead of failing the batch. Use to triage or rank a
    set of candidate structures.

    Args:
        params (BatchComputeMolecularPropertiesInput): smiles_list (1–100), response_format.

    Returns:
        Markdown summary table over all parseable molecules plus a list of any that failed.
    """
    results = []
    failed: list[str] = []
    for smi in params.smiles_list:
        props = _compute_molecular_properties(smi)
        if props is None:
            failed.append(smi)
        else:
            results.append(props)
    batch = BatchMolecularProperties(
        results=results, failed=failed, total_requested=len(params.smiles_list)
    )
    return _fmt(batch, params.response_format, "")


class CorpusDescribeInput(BaseModel):
    """Input for corpus_describe."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD


class CorpusSearchTargetsInput(BaseModel):
    """Input for corpus_search_targets_by_sequence."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    query: str = Field(
        ...,
        description="A target already in the corpus, given as an HGNC gene symbol (e.g. 'BRAF') "
        "or UniProt accession. The query is matched to its stored embedding — no sequence is "
        "embedded at request time.",
        min_length=1,
        max_length=50,
    )
    top_k: int = Field(
        default=10, ge=1, le=50, description="Maximum number of neighbor targets to return."
    )
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD


class CorpusFindSimilarCompoundsInput(BaseModel):
    """Input for corpus_find_similar_compounds."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    smiles: str = Field(
        ...,
        description="Query SMILES, e.g. 'CC(=O)Oc1ccccc1C(=O)O' (aspirin). The ECFP4 Morgan "
        "fingerprint is computed locally (RDKit) and matched against the corpus.",
        min_length=1,
        max_length=2000,
    )
    top_k: int = Field(
        default=20, ge=1, le=100, description="Maximum number of similar compounds to return."
    )
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD


class CorpusSearchCompoundsInput(BaseModel):
    """Input for corpus_search_compounds (hybrid relational + similarity search)."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")
    target: str | None = Field(
        None,
        description="Restrict to compounds active against this target — gene symbol or UniProt "
        "accession in the corpus (e.g. 'BRAF'). Omit to search across all targets.",
        max_length=50,
    )
    standard_type: Literal["IC50", "Ki", "Kd", "EC50"] | None = Field(
        None, description="Filter to a single assay endpoint."
    )
    min_pchembl: float | None = Field(
        None, ge=0, le=15, description="Minimum pChEMBL (−log10 molar potency); 6 ≈ 1 µM, 9 ≈ 1 nM."
    )
    max_mol_weight: float | None = Field(
        None, gt=0, description="Maximum molecular weight (Da), e.g. 500 for a lead-like filter."
    )
    similar_to_smiles: str | None = Field(
        None,
        description="Optional query SMILES — when given, results are ranked by ECFP4 Tanimoto "
        "to it (computed locally with RDKit); otherwise ranked by potency.",
        max_length=2000,
    )
    limit: int = Field(default=25, ge=1, le=100, description="Max compounds to return (page size).")
    offset: int = Field(default=0, ge=0, description="Pagination offset.")
    response_format: Literal["markdown", "json"] = _RESPONSE_FORMAT_FIELD


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=False
    )
)
async def corpus_describe(params: CorpusDescribeInput) -> str:
    """Describe the embedding-backed corpus store: coverage, counts, and provenance.

    The corpus is a curated, indexed bioactivity store (Postgres + pgvector) that the other
    corpus_* tools query — distinct from the live-API tools. Call this FIRST to establish what
    the corpus covers (target family, record counts) and from which upstream releases it was
    built (ChEMBL release, UniProt snapshot, embedding models) before trusting any corpus
    retrieval. Closed-world: results come only from the indexed corpus, not live APIs.

    Args:
        params (CorpusDescribeInput): response_format.

    Returns:
        Markdown manifest (target family, target/compound/activity counts, source versions,
        embedding models, build date). Reports an error if the corpus is not configured
        (no GENESIS_CORPUS_DSN) or has not been built yet.
    """
    pool = getattr(mcp.state, "corpus_pool", None)
    if pool is None:
        return _fmt(
            None,
            params.response_format,
            "Corpus store is not configured. Set GENESIS_CORPUS_DSN to a Postgres+pgvector "
            "instance and build the corpus (see docs/ROADMAP.md v0.6.0).",
            error_status="UpstreamUnavailable",
        )
    manifest = await fetch_manifest(pool)
    if manifest is None:
        return _fmt(
            None,
            params.response_format,
            "Corpus store is configured but has not been built yet — run the ingestion "
            "pipeline to populate it.",
            error_status="NotFound",
        )
    result = CorpusManifest(
        target_family=manifest["target_family"],
        built_at=manifest["built_at"].isoformat() if manifest.get("built_at") else "",
        chembl_release=manifest.get("chembl_release"),
        uniprot_snapshot=manifest.get("uniprot_snapshot"),
        protein_embedding_model=manifest.get("protein_embedding_model"),
        chem_embedding_model=manifest.get("chem_embedding_model"),
        target_count=manifest.get("target_count", 0),
        compound_count=manifest.get("compound_count", 0),
        activity_count=manifest.get("activity_count", 0),
    )
    return _fmt(result, params.response_format, "")


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=False
    )
)
async def corpus_search_targets_by_sequence(params: CorpusSearchTargetsInput) -> str:
    """Find indexed targets most similar to a query target by ESM-2 sequence embedding.

    Closed-world similarity search over the corpus: give a target already in the corpus (gene
    symbol or UniProt accession), and this returns the other corpus targets whose protein
    sequence embedding is closest (pgvector cosine kNN), each with a bioactivity-coverage count.
    Useful for finding sequence-related targets (e.g. kinases in the same subfamily) that share
    chemical matter. No sequence is embedded at request time — the query is matched to its
    stored vector, so this never loads an ML model. Call corpus_describe first to see coverage.

    Args:
        params (CorpusSearchTargetsInput): query (gene symbol / UniProt accession), top_k,
            response_format.

    Returns:
        Markdown table of neighbor targets ranked by cosine similarity (with kinase group and
        activity counts). Reports an error if the corpus is unconfigured, or if the query target
        is not in the corpus. Similarity is a retrieval signal, not a validated relationship.
    """
    pool = getattr(mcp.state, "corpus_pool", None)
    if pool is None:
        return _fmt(
            None,
            params.response_format,
            "Corpus store is not configured. Set GENESIS_CORPUS_DSN and build the corpus "
            "(see docs/ROADMAP.md v0.6.0).",
            error_status="UpstreamUnavailable",
        )
    found = await search_similar_targets(pool, params.query, params.top_k)
    if found is None:
        return _fmt(
            None,
            params.response_format,
            f"'{params.query}' is not in the indexed corpus (the human kinome) or has no "
            "embedding. Use corpus_describe to see what is covered.",
            error_status="NotFound",
        )
    query_gene, hits = found
    result = CorpusTargetNeighbors(
        query_gene=query_gene,
        hits=[CorpusTargetHit(**hit) for hit in hits],
    )
    return _fmt(result, params.response_format, "")


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=False
    )
)
async def corpus_find_similar_compounds(params: CorpusFindSimilarCompoundsInput) -> str:
    """Find indexed compounds most similar to a query molecule (ECFP4 Tanimoto).

    Closed-world chemical-similarity search over the corpus: the query SMILES is fingerprinted
    locally with RDKit (Morgan/ECFP4 — no ML model), then matched against the corpus by
    Tanimoto (pgvector Jaccard over `bit(2048)` fingerprints). Returns analog compounds with
    their Tanimoto and a count of their bioactivity records in the corpus — useful for finding
    known chemical matter near a hit, or analogs active against the indexed kinome. Call
    corpus_describe first to see what is covered.

    Args:
        params (CorpusFindSimilarCompoundsInput): smiles, top_k, response_format.

    Returns:
        Markdown table of similar compounds ranked by Tanimoto (with MW and activity counts).
        Reports an error if the corpus is unconfigured, the SMILES is invalid, or no
        fingerprinted compounds are indexed. Similarity is a retrieval signal, not a validated
        activity prediction.
    """
    pool = getattr(mcp.state, "corpus_pool", None)
    if pool is None:
        return _fmt(
            None,
            params.response_format,
            "Corpus store is not configured. Set GENESIS_CORPUS_DSN and build the corpus "
            "(see docs/ROADMAP.md v0.6.0).",
            error_status="UpstreamUnavailable",
        )
    fp_bits = _morgan_fp_bits(params.smiles)
    if fp_bits is None:
        return _fmt(
            None,
            params.response_format,
            f"Invalid SMILES: '{params.smiles}' could not be parsed by RDKit.",
            error_status="InvalidInput",
        )
    rows = await _corpus_search_similar_compounds(pool, fp_bits, params.top_k)
    if not rows:
        return _fmt(
            None,
            params.response_format,
            "No fingerprinted compounds are indexed yet — run the compound ingestion "
            "(see docs/ROADMAP.md v0.6.0).",
            error_status="NotFound",
        )
    result = CorpusSimilarCompounds(
        query_smiles=params.smiles,
        hits=[CorpusCompoundHit(**row) for row in rows],
    )
    return _fmt(result, params.response_format, "")


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=False
    )
)
async def corpus_search_compounds(params: CorpusSearchCompoundsInput) -> str:
    """Hybrid search of the corpus: relational filters + optional structural similarity, in one query.

    The headline corpus tool. Combine SQL-style filters — target (gene/UniProt), assay type,
    minimum pChEMBL potency, maximum molecular weight — with optional Morgan/Tanimoto similarity
    to a query SMILES, and get back bioactive compounds ranked by similarity (if a query molecule
    is given) or by potency. This answers questions neither a database nor a vector index can
    alone, e.g. "compounds structurally like this hit that ALSO have sub-µM measured activity
    against this kinase". Each hit carries three honestly-separate signals: Tanimoto (retrieval
    closeness), assay confidence (data quality), and pChEMBL (potency). Closed-world over the
    corpus; the query fingerprint is computed locally with RDKit (no ML model).

    Args:
        params (CorpusSearchCompoundsInput): target, standard_type, min_pchembl, max_mol_weight,
            similar_to_smiles, limit, offset, response_format.

    Returns:
        Markdown table of matching compounds with potency, assay confidence, MW (and Tanimoto in
        similarity mode). Errors if the corpus is unconfigured, the query SMILES is invalid, or
        the named target is not in the corpus. Similarity is a retrieval signal, not a validated
        activity prediction.
    """
    pool = getattr(mcp.state, "corpus_pool", None)
    if pool is None:
        return _fmt(
            None,
            params.response_format,
            "Corpus store is not configured. Set GENESIS_CORPUS_DSN and build the corpus "
            "(see docs/ROADMAP.md v0.6.0).",
            error_status="UpstreamUnavailable",
        )
    fp_bits: str | None = None
    if params.similar_to_smiles:
        fp_bits = _morgan_fp_bits(params.similar_to_smiles)
        if fp_bits is None:
            return _fmt(
                None,
                params.response_format,
                f"Invalid SMILES: '{params.similar_to_smiles}' could not be parsed by RDKit.",
                error_status="InvalidInput",
            )
    accession: str | None = None
    if params.target:
        accession = await resolve_target_accession(pool, params.target)
        if accession is None:
            return _fmt(
                None,
                params.response_format,
                f"Target '{params.target}' is not in the indexed corpus (the human kinome). "
                "Use corpus_describe to see coverage.",
                error_status="NotFound",
            )
    rows = await hybrid_search_compounds(
        pool,
        uniprot_accession=accession,
        standard_type=params.standard_type,
        min_pchembl=params.min_pchembl,
        max_mol_weight=params.max_mol_weight,
        query_fp_bits=fp_bits,
        limit=params.limit,
        offset=params.offset,
    )
    result = CorpusCompoundSearchResults(
        target=params.target,
        standard_type=params.standard_type,
        min_pchembl=params.min_pchembl,
        max_mol_weight=params.max_mol_weight,
        similar_to_smiles=params.similar_to_smiles,
        total=len(rows),
        hits=[CorpusCompoundSearchHit(**row) for row in rows],
    )
    return _fmt(result, params.response_format, "")


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
    if not __import__("os").environ.get("ANTHROPIC_API_KEY"):
        return (
            "**run_biology_workflow unavailable:** `ANTHROPIC_API_KEY` is not set.\n\n"
            "This tool drives an internal Claude agent to chain database queries, which "
            "requires an Anthropic API key. Set `ANTHROPIC_API_KEY=<your-key>` in the server "
            "environment (e.g. under `env` in claude_desktop_config.json), then restart the "
            "server. The individual tools (e.g. `prioritize_target`, `get_target_disease_association`) "
            "work without it."
        )
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
        ``use_when`` fields for every registered tool (count is sourced
        from ``len(registry)`` so it stays correct as tools are added).
    """
    registry = build_tool_registry(mcp.state)
    return format_registry_docs(registry)


@mcp.resource("health://status")
async def health_status_resource() -> str:
    """Live server health: upstream reachability, cache state, version, and uptime.

    An ops/observability surface (mirroring the served model's own ``/health``) that lets a
    client, orchestrator, or agent verify the server and its critical upstream services are
    reachable before relying on them. Probes a representative set of upstreams concurrently
    (each with a short timeout) and reports per-service reachability + latency.

    Returns:
        Markdown with overall status, version, tool count, uptime, DepMap cache size, and an
        upstream-reachability table.
    """
    state = mcp.state
    tools = await mcp.list_tools()
    uptime = time.time() - getattr(state, "started_at", time.time())
    client = getattr(state, "http_client", None)
    upstreams = await check_upstreams(client) if client is not None else []
    report = build_health_report(
        n_tools=len(tools),
        uptime_s=uptime,
        depmap_gene_count=getattr(state, "depmap_gene_count", 0),
        upstreams=upstreams,
    )
    return health_to_markdown(report)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
