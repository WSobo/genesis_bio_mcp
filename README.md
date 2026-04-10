# genesis-bio-mcp

[![CI](https://github.com/WSobo/genesis-bio-mcp/actions/workflows/ci.yml/badge.svg)](https://github.com/WSobo/genesis-bio-mcp/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.11%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![MCP](https://img.shields.io/badge/protocol-MCP-purple)](https://modelcontextprotocol.io)

An MCP (Model Context Protocol) server that connects AI agents to major public biomedical databases for drug discovery target prioritization and experiment design.

Ask an AI agent *"Assess BRAF as an oncology target for melanoma"* and watch it autonomously chain queries across UniProt, Open Targets, DepMap, GWAS Catalog, ChEMBL, AlphaFold, STRING, DGIdb, ClinicalTrials.gov, and Reactome into a structured evidence report — no hardcoded scripts, no manual API calls.

## Tools

### Core — Target Prioritization

| Tool | Database(s) | Purpose |
|------|-------------|---------|
| `resolve_gene` | UniProt + NCBI | Resolve gene aliases → canonical HGNC symbol, NCBI ID, UniProt accession |
| `get_protein_info` | UniProt Swiss-Prot | Protein function, pathways, disease variants, subcellular location |
| `get_target_disease_association` | Open Targets | Evidence-based association score (0–1) for a target–disease pair |
| `get_cancer_dependency` | DepMap + Open Targets | CRISPR essentiality scores across cancer cell lines |
| `get_gwas_evidence` | GWAS Catalog | Genome-wide significant SNP associations for a trait |
| `get_compounds` | PubChem + ChEMBL | Active small molecules with potency data (IC50/Ki/Kd) |
| `prioritize_target` | All of the above | Full parallel evidence synthesis → priority score (0–10) + report |
| `compare_targets` | All of the above | Rank 2–5 targets side by side for an indication |

### Lab Loop — Experiment Design

| Tool | Database(s) | Purpose |
|------|-------------|---------|
| `get_protein_structure` | AlphaFold + RCSB PDB | Structure confidence, resolution, co-crystallized ligands |
| `get_protein_interactome` | STRING | High-confidence interaction partners (selectivity risk profile) |
| `get_drug_history` | DGIdb + ClinicalTrials.gov | Known drugs, approval status, trial counts by phase |
| `get_pathway_context` | Reactome | Enriched pathway membership with statistical significance |

## Quickstart

### Install

```bash
git clone https://github.com/WSobo/genesis-bio-mcp
cd genesis-bio-mcp
uv sync
```

### Add to Claude Desktop

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

### Try it

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
├── server.py                      # FastMCP server, tool registration, lifespan
├── models.py                      # Pydantic V2 output models (all with to_markdown())
├── clients/
│   ├── uniprot.py                 # UniProt REST: Swiss-Prot gene + protein data
│   ├── open_targets.py            # Open Targets GraphQL: 3-step target–disease resolution
│   ├── depmap.py                  # DepMap task API + disk cache + OT lineage fallback
│   ├── gwas.py                    # GWAS Catalog HAL/REST, Unicode normalization
│   ├── pubchem.py                 # PubChem REST, Semaphore rate limiting, tenacity retries
│   ├── chembl.py                  # ChEMBL REST: target lookup + pChEMBL potency data
│   ├── alphafold.py               # AlphaFold + RCSB PDB structural data
│   ├── string_db.py               # STRING protein interaction network
│   ├── dgidb.py                   # DGIdb 5.0 GraphQL: known drugs targeting gene
│   ├── clinical_trials.py         # ClinicalTrials.gov v2 REST
│   └── reactome.py                # Reactome pathway enrichment (two-step token API)
└── tools/
    ├── gene_resolver.py           # Multi-source alias resolution
    └── target_prioritization.py   # asyncio.gather orchestration, score computation
```

### Key design decisions

| Decision | Rationale |
|----------|-----------|
| Single shared `httpx.AsyncClient` via lifespan | Connection pooling across all 12 tools |
| `asyncio.gather` for all sub-queries | All APIs hit simultaneously; no sequential blocking |
| `_safe()` wrapping on every coroutine | One API failure never crashes a pipeline |
| Pydantic V2 models with `to_markdown()` | Tools output agent-readable markdown, never raw JSON |
| DepMap disk cache (7-day TTL) | Warm starts in <1s vs 30–60s cold download |
| Extended mode via `extended=True` | Lab loop tools run in the same gather; score unchanged |
| `asyncio.Semaphore` per rate-limited API | STRING/ChEMBL: `Semaphore(2)`, Reactome: `Semaphore(3)` |

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
3. Register a `@mcp.tool()` in `server.py` — the tool must return a markdown string, not a dict
4. Inject the client via `lifespan` and pass through `server.state`
5. Add unit tests in `tests/test_clients.py` using `respx` to mock HTTP at the transport level
6. Run `uv run ruff format . && uv run ruff check --fix .` and `uv run pytest tests/ -v`
7. Update this README

## Environment variables

| Variable | Default | Description |
|----------|---------|-------------|
| `NCBI_EMAIL` | `genesis-bio-mcp@example.com` | Required for NCBI E-utils polite use |

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

## License

MIT
