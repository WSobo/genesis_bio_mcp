# genesis-bio-mcp

[![CI](https://github.com/WSobo/genesis-bio-mcp/actions/workflows/ci.yml/badge.svg)](https://github.com/WSobo/genesis-bio-mcp/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.11%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![MCP](https://img.shields.io/badge/protocol-MCP-purple)](https://modelcontextprotocol.io)

An MCP server that connects AI agents to 12 biomedical databases for drug discovery target prioritization and experiment design.

Ask *"Find underexplored MAPK kinases with no approved drugs"* and a Claude-powered workflow agent chains queries across UniProt, Open Targets, DepMap, GWAS Catalog, ChEMBL, PubChem, AlphaFold, STRING, DGIdb, ClinicalTrials.gov, and Reactome into a structured evidence report — no hardcoded scripts, no manual API calls.

> **What is MCP?** The [Model Context Protocol](https://modelcontextprotocol.io) is an open standard that lets AI assistants (Claude, Cursor, etc.) call external tools and data sources. This server exposes biomedical databases as MCP tools so any MCP-compatible AI can query them directly during a conversation.

## Install

**uv is required** ([installation guide](https://docs.astral.sh/uv/getting-started/installation/)):

```bash
# macOS / Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# Windows
winget install --id=astral-sh.uv
```

Then clone and sync:

```bash
git clone https://github.com/WSobo/genesis-bio-mcp
cd genesis-bio-mcp
uv sync
```

> pip and conda are not supported — the project uses uv for reproducible dependency resolution.

For the AI workflow tool (`run_biology_workflow`) only:

```bash
export ANTHROPIC_API_KEY=sk-ant-...
```

## Add to Claude Desktop

Edit `~/Library/Application Support/Claude/claude_desktop_config.json` (macOS) or `%APPDATA%\Claude\claude_desktop_config.json` (Windows):

```json
{
  "mcpServers": {
    "genesis-bio-mcp": {
      "command": "uv",
      "args": ["run", "--directory", "/path/to/genesis-bio-mcp", "genesis-bio-mcp"]
    }
  }
}
```

Replace `/path/to/genesis-bio-mcp` with the absolute path to your clone. Restart Claude Desktop.

## Tools

### Core — Target Prioritization

| Tool | Database(s) | Purpose |
|------|-------------|---------|
| `resolve_gene` | UniProt + NCBI | Resolve gene aliases → canonical HGNC symbol, NCBI ID, UniProt accession |
| `get_protein_info` | UniProt Swiss-Prot | Protein function, pathways, disease variants, subcellular location |
| `get_target_disease_association` | Open Targets | Evidence-based association score (0–1) for a target–disease pair |
| `get_cancer_dependency` | DepMap | CRISPR essentiality scores across cancer cell lines |
| `get_gwas_evidence` | GWAS Catalog | Genome-wide significant SNP associations for a trait |
| `get_compounds` | PubChem | Active small molecules with bioactivity hit counts |
| `get_chembl_compounds` | ChEMBL | Quantitative potency data (IC50/Ki/Kd) with pChEMBL values |
| `prioritize_target` | All of the above | Full parallel evidence synthesis → priority score (0–10) + confidence interval |
| `compare_targets` | All of the above | Rank 2–5 targets side by side for an indication |

### Lab Loop — Experiment Design

| Tool | Database(s) | Purpose |
|------|-------------|---------|
| `get_protein_structure` | AlphaFold + RCSB PDB | Structure confidence (pLDDT), resolution, co-crystallized ligands |
| `get_protein_interactome` | STRING | High-confidence interaction partners (selectivity risk profile) |
| `get_drug_history` | DGIdb + ClinicalTrials.gov | Known drugs, approval status, trial counts by phase |
| `get_pathway_context` | Reactome | Enriched pathway membership with statistical significance |

### AI Workflow — Dynamic Orchestration

| Tool | Requires | Purpose |
|------|----------|---------|
| `run_biology_workflow` | `ANTHROPIC_API_KEY` | Free-text question → Claude selects and chains tools autonomously |

Unlike the other tools, `run_biology_workflow` does not have a fixed pipeline. It instantiates an inner Claude agent (`claude-sonnet-4-6`) that reads your question, reasons about which databases are relevant, calls tools (in parallel when possible), and synthesizes a narrative answer with citations to the data it retrieved. Useful for open-ended questions where the right combination of databases isn't known in advance.

All tools accept a `response_format` parameter (`"markdown"` default, `"json"` for programmatic integration).

## Try it

**Target prioritization:**
- *"Is PCSK9 a good target for cardiovascular disease? Check Open Targets and PubChem."*
- *"Assess BRAF as an oncology target for melanoma — full report."*
- *"What GWAS evidence links FTO to obesity?"*
- *"Compare BRAF, EGFR, and KRAS for non-small cell lung cancer."*

**Lab loop / experiment design:**
- *"Does BRAF have a co-crystal structure with an inhibitor? What's the resolution?"*
- *"What proteins does EGFR interact with? Are any in the same family?"*
- *"What drugs already target KRAS? Are any approved?"*
- *"What signaling pathways is PIK3CA part of? Any combination opportunities?"*

**AI workflow (requires `ANTHROPIC_API_KEY`):**
- *"Find underexplored kinases in the MAPK pathway with no approved drugs."*
- *"Compare BRAF and MEK1 for melanoma — which has better chemical matter?"*
- *"What is the target validation evidence for FTO in obesity?"*

## Example Output

### Standard mode

```
prioritize_target("BRAF", "melanoma")

→ priority_score: 7.11 / 10
→ priority_tier: High
→ evidence_summary: "BRAF shows strong Open Targets association with melanoma (score: 0.82,
  n=5 evidence items). Open Targets reports strong known-drug evidence (score: 0.98),
  suggesting existing approved therapeutics. DepMap CRISPR data show dependency in 9% of
  cancer lines, highest in cancer or benign tumor. ChEMBL reports 68 compounds; best IC50 ≈
  0.3 nM (pChEMBL=9.5)."

→ disease_association.overall_score:          0.82
→ disease_association.somatic_mutation_score: 0.80   # BRAF V600E is a major somatic driver
→ disease_association.known_drug_score:       0.98   # vemurafenib, dabrafenib, encorafenib
→ cancer_dependency.fraction_dependent_lines: 0.095  # 115/1208 lines — selective
→ chembl_compounds.best_pchembl:              9.52
→ chembl_compounds.total_active_compounds:    68
→ data_gaps: ["gwas"]
```

**Score breakdown:** OT (0.82×3=2.46) + DepMap (0.095×2=0.19) + clinical (0.98×1.5=1.46) + ChEMBL pChEMBL≥9 (1.5) + protein (reviewed+variants=1.5) = **7.11**.

> **Note on GWAS for BRAF/melanoma**: BRAF V600E is a somatic driver, not a germline variant. GWAS Catalog correctly returns no melanoma hits near BRAF. Reported as `data_gap`. GWAS is strong for germline-driven targets like *FTO* (obesity) or *PCSK9* (cardiovascular disease).

**Confidence Assessment:**

```
→ data_coverage_pct:         83.3   # 5 of 6 core sources returned data
→ score_confidence_interval: (6.0, 8.2)  # ±0.7 due to missing GWAS
→ proxy_data_flags:          {}     # all data is real (no OT proxies used)
```

The CI widens automatically when data sources are unavailable. For BRAF/melanoma the missing GWAS is expected (somatic driver) — the CI correctly signals reduced certainty rather than silently accepting the gap.

### Extended mode

Pass `extended=True` to include all four lab loop tools in a single gathered call:

```python
prioritize_target("BRAF", "melanoma", extended=True)
```

Adds to the report:

```
→ protein_structure.alphafold_plddt:    92.1       # high confidence (≥90)
→ protein_structure.best_resolution:    1.7 Å
→ protein_structure.has_ligand_bound:   true       # inhibitor co-crystal available
→ protein_structure.total_pdb:         156

→ protein_interactome.top_partners:    MAP2K1 (0.999), MAP2K2 (0.998), RAF1 (0.963)

→ drug_history.approved_drug_count:    4           # vemurafenib, dabrafenib, encorafenib, ...
→ drug_history.trial_counts_by_phase:  {"Phase 1": 12, "Phase 2": 8, "Phase 3": 3}

→ pathway_context.top_pathway:         "MAPK1/MAPK2 Cascade" (p=2.3e-15)
→ pathway_context.category:            Signaling
```

## Scoring Model

The composite priority score (0–10) combines six evidence axes:

| Source | Max | Logic |
|--------|-----|-------|
| Open Targets association | 3.0 | `overall_score × 3` |
| DepMap CRISPR dependency | 2.0 | `fraction_dependent × 2` (×0.7 if OT proxy used) |
| GWAS evidence | 2.0 | `min(hits, 10) / 10 × 2` |
| Clinical / known-drug evidence | 1.5 | `known_drug_score × 1.5` |
| ChEMBL potency | 1.5 | pChEMBL ≥9 → 1.5, ≥7 → 1.0, ≥5 → 0.5, else 0.25 |
| UniProt protein quality | 1.5 | reviewed (+0.5) + variant coverage (max +1.0) |

Pan-essential genes (DepMap `common_essential`) have their DepMap contribution capped at 0.5.

## Architecture

```
src/genesis_bio_mcp/
├── server.py                      # FastMCP server, 14 tools, lifespan, tool://registry resource
├── models.py                      # Pydantic V2 output models (all with to_markdown() + model_dump_json())
├── workflow_agent.py              # ToolSpec registry, run_agent_loop(), format_registry_docs()
├── clients/
│   ├── uniprot.py                 # UniProt REST (session-scoped cache)
│   ├── open_targets.py            # Open Targets GraphQL: 3-step target–disease resolution
│   ├── depmap.py                  # DepMap task API + disk cache + OT lineage fallback
│   ├── gwas.py                    # GWAS Catalog HAL/REST, Unicode normalization
│   ├── pubchem.py                 # PubChem REST, Semaphore rate limiting, tenacity retries
│   ├── chembl.py                  # ChEMBL REST: target lookup + pChEMBL potency data
│   ├── alphafold.py               # AlphaFold + RCSB PDB (session-scoped cache)
│   ├── string_db.py               # STRING protein interaction network
│   ├── dgidb.py                   # DGIdb 5.0 GraphQL
│   ├── clinical_trials.py         # ClinicalTrials.gov v2 REST
│   └── reactome.py                # Reactome pathway enrichment, two-step token API (session-scoped cache)
└── tools/
    ├── gene_resolver.py           # Multi-source alias resolution (HER2→ERBB2, p53→TP53)
    └── target_prioritization.py   # asyncio.gather orchestration, score computation, confidence CI
```

### Key design decisions

| Decision | Rationale |
|----------|-----------|
| Single shared `httpx.AsyncClient` via lifespan | Connection pooling across all 14 tools |
| `asyncio.gather` for all sub-queries | All APIs hit simultaneously; no sequential blocking |
| `_safe()` wrapping on every coroutine | One API failure never crashes a pipeline |
| Pydantic V2 models with `to_markdown()` | Tools output agent-readable markdown, never raw JSON |
| DepMap disk cache (7-day TTL) | Warm starts in <1s vs 30–60s cold download |
| Extended mode via `extended=True` | Lab loop tools run in the same gather; score unchanged |
| `asyncio.Semaphore` per rate-limited API | STRING/ChEMBL: `Semaphore(2)`, Reactome: `Semaphore(3)` |
| `workflow_agent.py` with `ToolSpec` registry | Each tool has `tool_category` + embedding-searchable `use_when` — enables dynamic selection without hardcoded routing |
| Session-scoped dict cache on UniProt, AlphaFold, Reactome | Repeated queries (e.g. `compare_targets` on the same gene) served from memory; no redundant HTTP round-trips |
| `response_format` param on all tools | Default `"markdown"` for agent consumption; `"json"` for programmatic / pipeline integration |

### MCP Resource: `tool://registry`

Clients and agent frameworks can read `tool://registry` to get a structured Markdown catalogue of all 14 tools, grouped by category, with `use_when` guidance for semantic retrieval. No tool call needed — useful for embedding-based tool selection or for auditing available capabilities.

## Development

```bash
uv sync
uv run pytest tests/ -v
uv run pytest tests/ -v --cov=genesis_bio_mcp

# Full integration test against live APIs (saves reports to examples/)
uv run python test_full.py

# Single target
uv run python test_full.py BRAF melanoma

# Extended mode (lab loop tools)
uv run python test_full.py BRAF melanoma --extended
uv run python test_full.py PCSK9 hypercholesterolemia --extended
```

### Running the server

```bash
uv run genesis-bio-mcp          # stdio transport (MCP clients)
```

### Lint & format

```bash
uv run ruff format .
uv run ruff check --fix .
```

### Adding a new data source

1. Create `src/genesis_bio_mcp/clients/<source>.py` — wrap HTTP calls in try/except, return `None` on error
2. Add Pydantic V2 output model to `models.py` with `to_markdown()` returning a markdown string
3. Register a `@mcp.tool()` in `server.py` — the tool must return a markdown string, not a dict (14 tools currently)
4. Inject the client via `lifespan` and pass through `server.state`
5. Add a `ToolSpec` entry in `workflow_agent.py`'s `build_tool_registry()` with `tool_category` and `use_when`
6. Add unit tests in `tests/test_clients.py` using `respx` to mock HTTP at the transport level
7. Run `uv run ruff format . && uv run ruff check --fix .` and `uv run pytest tests/ -v`
8. Update this README

## Environment variables

| Variable | Default | Required for | Description |
|----------|---------|--------------|-------------|
| `ANTHROPIC_API_KEY` | — | `run_biology_workflow` | Anthropic API key for the inner Claude agent |
| `NCBI_EMAIL` | `genesis-bio-mcp@example.com` | All NCBI queries | Required for NCBI E-utils polite use |

## API reference

| Database | API Type | Rate Limit | Notes |
|----------|----------|------------|-------|
| UniProt | REST | Generous | `organism_id:9606 AND reviewed:true` for human Swiss-Prot |
| Open Targets | GraphQL | ~2 req/s | Requires Ensembl ID + EFO ID — resolved automatically |
| DepMap | REST (task queue) | Moderate | Celery polling → pre-signed CSV URL. 7-day disk cache. OT fallback if unreachable |
| GWAS Catalog | REST/HAL | Moderate | Unicode normalization required for trait matching |
| PubChem | REST | 5 req/s | HTTP 503 on rate limit; tenacity + Semaphore |
| ChEMBL | REST | ~1 req/s | No API key. Two-step: target → bioactivity (IC50/Ki/Kd). `Semaphore(2)` |
| AlphaFold | REST | Generous | pLDDT ≥90 = high confidence; <70 = disordered |
| RCSB PDB | REST | Generous | Search → entry fetch per structure |
| STRING | REST | ~2 req/s | `required_score=700` for high-confidence interactions. `Semaphore(2)` |
| DGIdb | GraphQL | Generous | DGIdb 5.0; deduplicates drugs, surfaces approval status |
| ClinicalTrials.gov | REST v2 | Generous | Phase normalization: `PHASE1` → `"Phase 1"` |
| Reactome | REST | ~5 req/s | Two-step token API. `Semaphore(3)` |
| NCBI E-utils | REST | 3 req/s | Set `NCBI_EMAIL` env var |

