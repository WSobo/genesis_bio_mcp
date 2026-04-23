# genesis_bio_mcp

[![CI](https://github.com/WSobo/genesis-bio-mcp/actions/workflows/ci.yml/badge.svg)](https://github.com/WSobo/genesis-bio-mcp/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.11%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![MCP](https://img.shields.io/badge/protocol-MCP-purple)](https://modelcontextprotocol.io)

An MCP server that gives AI agents structured access to **23 tools** across
**19 biomedical databases** for drug discovery target prioritization,
experiment design, and protein engineering.

Ask *"Find underexplored MAPK kinases with no approved drugs"* or
*"Is TP53 R175H pathogenic — what does AlphaMissense and DMS say?"* and a
Claude-powered workflow agent chains queries across UniProt, Open Targets,
DepMap, GWAS Catalog, ChEMBL, PubChem, AlphaFold, STRING, DGIdb,
ClinicalTrials.gov, Reactome, SAbDab, IEDB, InterPro, MaveDB, gnomAD, and
MyVariant.info into a structured evidence report. no hardcoded scripts,
no manual API calls.

> **What is MCP?** The [Model Context Protocol](https://modelcontextprotocol.io)
> is an open standard that lets AI assistants (Claude, Cursor, etc.) call
> external tools and data sources. This server exposes biomedical
> databases as MCP tools so any MCP-compatible AI can query them directly
> during a conversation.

---

## Documentation

| Doc | What it covers |
|---|---|
| **[docs/tools.md](docs/tools.md)** | Full catalog of all 23 tools grouped by category, input fields, and use cases |
| **[docs/protein-engineering.md](docs/protein-engineering.md)** | v0.2.0 protein engineering workflows: sequence analysis, variant effects, T-cell immunogenicity, combined examples |
| **[docs/architecture.md](docs/architecture.md)** | Directory layout, client/model patterns, design decisions, `prioritize_target` scoring model, per-database API reference |
| **[docs/benchmark.md](docs/benchmark.md)** | 12-target benchmark matrix and example `prioritize_target` output |
| **[docs/deployment.md](docs/deployment.md)** | Environment variables, Docker setup, Claude Desktop config, production checklist |
| **[CONTRIBUTING.md](CONTRIBUTING.md)** | 7-step walkthrough for adding a new data source |
| **[CHANGELOG.md](CHANGELOG.md)** | Version history |

---

## Install

**uv is required** ([install guide](https://docs.astral.sh/uv/getting-started/installation/)):

```bash
# macOS / Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# Windows
winget install --id=astral-sh.uv
```

```bash
git clone https://github.com/WSobo/genesis-bio-mcp
cd genesis-bio-mcp
uv sync
```

> pip and conda are not supported — the project uses uv for reproducible
> dependency resolution.

For the AI workflow tool (`run_biology_workflow`) only:

```bash
export ANTHROPIC_API_KEY=sk-ant-...
```

## Add to Claude Desktop

Edit `~/Library/Application Support/Claude/claude_desktop_config.json` (macOS)
or `%APPDATA%\Claude\claude_desktop_config.json` (Windows):

```json
{
  "mcpServers": {
    "genesis-bio-mcp": {
      "command": "uv",
      "args": ["run", "--directory", "/path/to/genesis-bio-mcp", "genesis-bio-mcp"],
      "env": {
        "ANTHROPIC_API_KEY": "sk-ant-..."
      }
    }
  }
}
```

Replace `/path/to/genesis-bio-mcp` with the absolute path to your clone and
restart Claude Desktop. See [docs/deployment.md](docs/deployment.md) for
the full list of environment variables, Docker setup, and production
configuration.

---

## Tools at a glance

23 tools across 19 data sources, organized into 8 categories. Full details
in [docs/tools.md](docs/tools.md).

| Category | Tools | Data sources |
|---|---|---|
| **Gene annotation** | `resolve_gene`, `get_protein_info`, `get_protein_sequence` | UniProt, NCBI |
| **Disease evidence** | `get_target_disease_association`, `get_cancer_dependency`, `get_gwas_evidence` | Open Targets, DepMap, GWAS Catalog, EFO |
| **Druggability** | `get_compounds`, `get_chembl_compounds` | PubChem, ChEMBL |
| **Structure & interactions** | `get_protein_structure`, `get_protein_interactome`, `get_biogrid_interactions` | AlphaFold, RCSB PDB, STRING, BioGRID |
| **Antibody & epitope** | `get_antibody_structures`, `get_epitope_data`, `get_mhc_binding` | SAbDab, IEDB, IEDB NextGen Tools |
| **Protein engineering** | `get_protein_sequence`, `get_variant_effects`, `get_variant_constraints`, `get_domain_annotation`, `get_dms_scores`, `get_mhc_binding` | UniProt, gnomAD, MyVariant.info, MaveDB, InterPro, IEDB |
| **Pathways** | `get_pathway_context`, `get_pathway_members` | Reactome |
| **Synthesis** | `get_drug_history`, `prioritize_target`, `compare_targets`, `run_biology_workflow` | DGIdb + ClinicalTrials.gov + Claude |

All tools return Markdown by default; every tool accepts
`response_format="json"` for pipeline integration.

---

## Try it

### Target prioritization
- *"Is PCSK9 a good target for cardiovascular disease? Check Open Targets and PubChem."*
- *"Assess BRAF as an oncology target for melanoma — full report."*
- *"Compare BRAF, EGFR, and KRAS for non-small cell lung cancer."*

### Experiment design
- *"Does BRAF have a co-crystal structure with an inhibitor?"*
- *"What proteins does EGFR interact with?"*
- *"What drugs already target KRAS, and what stage are they in?"*

### Protein engineering
- *"What is the pathogenicity of TP53 R175H? Include AlphaMissense and DMS."*
- *"Get the BRAF kinase-domain sequence (residues 450–620) and flag deamidation hotspots."*
- *"Predict HLA-I T-cell epitopes in the peptide SLYNTVATL against the default human panel."*

### AI workflow (requires `ANTHROPIC_API_KEY`)
- *"Find underexplored kinases in the MAPK pathway with no approved drugs."*
- *"Is variant X in gene Y likely to break function? Check all available pathogenicity and DMS data."*
- *"Screen the CDR3 sequence ARDYRLDY for T-cell immunogenicity risk across common HLA alleles."*

---

## Development

```bash
uv sync
uv run pytest tests/ -v          # 146 unit + integration tests
uv run pytest tests/ --cov=genesis_bio_mcp
```

```bash
uv run genesis-bio-mcp           # stdio transport for MCP clients

uv run ruff format .
uv run ruff check --fix .
```

Full dev workflow and the 7-step recipe for adding a new data source:
[CONTRIBUTING.md](CONTRIBUTING.md).

---

## License

MIT
