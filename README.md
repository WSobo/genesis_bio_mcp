# genesis-bio-mcp

[![CI](https://github.com/WSobo/genesis-bio-mcp/actions/workflows/ci.yml/badge.svg)](https://github.com/WSobo/genesis-bio-mcp/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.11%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![MCP](https://img.shields.io/badge/protocol-MCP-purple)](https://modelcontextprotocol.io)

An MCP server that gives AI agents structured access to 14 tools across 12 biomedical databases for drug discovery target prioritization and experiment design.

Ask *"Find underexplored MAPK kinases with no approved drugs"* and a Claude-powered workflow agent chains queries across UniProt, Open Targets, DepMap, GWAS Catalog, ChEMBL, PubChem, AlphaFold, STRING, DGIdb, ClinicalTrials.gov, and Reactome into a structured evidence report — no hardcoded scripts, no manual API calls.

> **What is MCP?** The [Model Context Protocol](https://modelcontextprotocol.io) is an open standard that lets AI assistants (Claude, Cursor, etc.) call external tools and data sources. This server exposes biomedical databases as MCP tools so any MCP-compatible AI can query them directly during a conversation.

---

## Install

**uv is required** ([installation guide](https://docs.astral.sh/uv/getting-started/installation/)):

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

---

## Tools

### Core — Target Prioritization

| Tool | Database(s) | Purpose |
|------|-------------|---------|
| `resolve_gene` | UniProt + NCBI | Resolve gene aliases → canonical HGNC symbol, NCBI ID, UniProt accession |
| `get_protein_info` | UniProt Swiss-Prot | Protein function, pathways, disease variants, subcellular location |
| `get_target_disease_association` | Open Targets | Evidence-based association score (0–1) for a target–disease pair |
| `get_cancer_dependency` | DepMap | CRISPR essentiality scores across cancer cell lines |
| `get_gwas_evidence` | GWAS Catalog + EFO | Genome-wide significant SNP associations; trait matching via EFO ontology |
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

Unlike the other tools, `run_biology_workflow` does not have a fixed pipeline. It instantiates an inner Claude agent (`claude-sonnet-4-6`) that reads your question, reasons about which databases are relevant, calls tools in parallel where possible, and synthesizes a narrative answer with citations to the data it retrieved.

All tools accept a `response_format` parameter (`"markdown"` default, `"json"` for programmatic integration).

---

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

---

## Example Output

### Standard mode

```
prioritize_target("BRAF", "melanoma")

→ priority_score: 7.2 / 10
→ priority_tier: High
→ evidence_summary: "BRAF shows strong Open Targets association with melanoma (score: 0.82,
  n=5 evidence items). Open Targets reports strong known-drug evidence (score: 0.98),
  suggesting existing approved therapeutics. DepMap CRISPR data show dependency in 9%
  of cancer lines, highest in differentiated thyroid carcinoma, glioblastoma multiforme,
  lung adenocarcinoma. ChEMBL reports 68 compounds; best IC50 ≈ 0.3 nM (pChEMBL=9.5)."

→ disease_association.overall_score:          0.82
→ disease_association.somatic_mutation_score: 0.80   # BRAF V600E is the canonical somatic driver
→ disease_association.known_drug_score:       0.98   # vemurafenib, dabrafenib, encorafenib
→ cancer_dependency.fraction_dependent_lines: 0.09   # 9% — lineage-selective (melanoma)
→ chembl_compounds.best_pchembl:              9.5    # clinical-grade potency
→ chembl_compounds.total_active_compounds:    68
→ data_gaps: ["gwas"]                                # expected — BRAF V600E is somatic, not germline
```

**Confidence assessment:**

```
→ data_coverage_pct:         83.3   # 5 of 6 core sources returned data
→ score_confidence_interval: (6.0, 8.2)
→ proxy_data_flags:          {}     # all real data, no OT proxies used
```

### Extended mode

Pass `extended=True` to include all four lab loop tools in the same parallel gather:

```python
prioritize_target("BRAF", "melanoma", extended=True)
```

```
→ protein_structure.alphafold_plddt:    92.1       # high confidence (≥90)
→ protein_structure.best_resolution:    1.7 Å
→ protein_structure.has_ligand_bound:   true       # inhibitor co-crystal available
→ protein_interactome.top_partners:    MAP2K1 (0.999), MAP2K2 (0.998), RAF1 (0.963)
→ drug_history.approved_drug_count:    4
→ drug_history.trial_counts_by_phase:  {"Phase 1": 12, "Phase 2": 8, "Phase 3": 3}
→ pathway_context.top_pathway:         "MAPK1/MAPK2 Cascade" (p=2.3e-15)
```

---

## Benchmark

End-to-end results against live APIs — 12 representative targets spanning oncology, metabolic disease, autoimmune, and pain biology (single session, warm DepMap/EFO cache):

| Gene input | Resolved | Indication | Score | Tier | GWAS | Time |
|---|---|---|---|---|---|---|
| BRAF | BRAF | melanoma | 7.2 | High | — somatic driver | 4.3s |
| EGFR | EGFR | non-small cell lung carcinoma | 7.5 | High | — | 2.2s |
| KRAS | KRAS | pancreatic cancer | 4.9 | Medium | — | 2.2s |
| **HER2** | **ERBB2** | breast cancer | 7.6 | High | — | 2.0s |
| PCSK9 | PCSK9 | hypercholesterolemia | **9.1** | High | ✓ 4 hits, p=4×10⁻²⁰ | 2.3s |
| FTO | FTO | obesity | 4.3 | Medium | timeout† | 14.4s |
| TNF | TNF | rheumatoid arthritis | **7.8** | High | ✓ 1 hit, p=9×10⁻²⁵ | 9.1s |
| PTGS2 | PTGS2 | inflammation | 7.7 | High | — | 7.2s |
| TP53 | TP53 | squamous cell carcinoma | 3.6 | Low | — | 3.8s |
| CD274 | CD274 | melanoma | 4.1 | Medium | — | 6.3s |
| **p53** | **TP53** | lung cancer | 3.4 | Low | — session cache | 3.3s |
| **COX2** | **PTGS2** | pain | **7.8** | High | — session cache | **1.9s** |

Bold gene inputs indicate alias resolution. GWAS gaps for BRAF/EGFR/KRAS are biologically expected — somatic cancer drivers are not GWAS loci. COX2/pain at 1.9s demonstrates the session gene cache: PTGS2 associations were already fetched for the inflammation query and are reused instantly.

†FTO has extensive GWAS signal (strong obesity associations with p<10⁻¹⁰⁰) but the Catalog API is slow for high-association genes; disk cache rescues repeat queries.

---

## Scoring Model

The composite priority score (0–10) combines six evidence axes:

| Source | Max | Logic |
|--------|-----|-------|
| Open Targets association | 3.0† | `overall_score × 3` |
| DepMap CRISPR dependency | 2.0 | `fraction_dependent × 2` (×0.7 if OT proxy; ×1.2 if indication matches top lineage) |
| GWAS evidence | 2.0 | `min(hits, 3) / 3 × 2` — saturates at 3 replicated hits |
| Clinical / known-drug evidence | 1.5 | `known_drug_score × 1.5` |
| ChEMBL potency | 1.5 | pChEMBL ≥9 → 1.5, ≥7 → 1.0, ≥5 → 0.5, else 0.25 |
| UniProt protein quality | 1.5 | reviewed (+0.5) + variant coverage (max +1.0) |

Pan-essential genes (DepMap `common_essential`) have their DepMap contribution capped at 0.5 — a narrow therapeutic window is a liability, not an asset.

†**Biologics floor**: when `known_drug_score > 0.9` AND both `genetic_association_score` and `somatic_mutation_score` are null, the OT contribution is floored at 3.25. This prevents pharmacologically-validated biologics targets (HER2/breast cancer, TNF/RA) from being underscored because OT's multi-datatype composite is depressed by evidence classes irrelevant to their mechanism. Targets with populated genetic or somatic scores (BRAF/melanoma) are unaffected.

Scores are documented in `target_prioritization.py` with explicit reasoning for each constant value — calibrated against the benchmark matrix and intentionally conservative to avoid overfitting to any single target class.

---

## Architecture

```
src/genesis_bio_mcp/
├── server.py                      # FastMCP server, 14 tools, lifespan, tool://registry resource
├── models.py                      # Pydantic V2 output models (all with to_markdown() + model_dump_json())
├── workflow_agent.py              # ToolSpec registry, run_agent_loop(), format_registry_docs()
├── config/
│   ├── efo_resolver.py            # OLS4 EFO ontology client: free-text query → EFO terms + synonyms
│   └── trait_synonyms.py          # GWAS trait matching: EFO URI → synonym expansion → fallback dict
├── clients/
│   ├── uniprot.py                 # UniProt REST (session-scoped in-memory cache)
│   ├── open_targets.py            # Open Targets GraphQL: 3-step target–disease resolution; 5xx retry
│   ├── depmap.py                  # DepMap task API + 7-day disk cache + OT lineage proxy fallback
│   ├── gwas.py                    # GWAS Catalog HAL/REST: concurrent fetch paths, session gene cache
│   ├── pubchem.py                 # PubChem REST, Semaphore rate limiting, tenacity retries
│   ├── chembl.py                  # ChEMBL REST: target lookup + pChEMBL potency data
│   ├── alphafold.py               # AlphaFold + RCSB PDB (session-scoped cache)
│   ├── string_db.py               # STRING protein interaction network
│   ├── dgidb.py                   # DGIdb 5.0 GraphQL
│   ├── clinical_trials.py         # ClinicalTrials.gov v2 REST
│   └── reactome.py                # Reactome pathway enrichment (session-scoped cache)
└── tools/
    ├── gene_resolver.py           # Multi-source alias resolution (HER2→ERBB2, p53→TP53)
    └── target_prioritization.py   # asyncio.gather orchestration, score computation, confidence CI
```

### Key design decisions

| Decision | Rationale |
|----------|-----------|
| Single shared `httpx.AsyncClient` via lifespan | Connection pooling across all 14 tools; no per-request TLS handshake overhead |
| `asyncio.gather` for all sub-queries | All 6 core APIs hit simultaneously — total latency ≈ slowest single call, not sum |
| `_safe()` wrapping on every coroutine | One API failure never crashes the pipeline; errors surface in `data_gaps`/`errors` |
| EFO ontology-backed GWAS trait matching | Free-text queries ("fat", "joint inflammation", "sugar disease") resolve via OLS4 to EFO terms. Matching uses EFO URI exact match (gene-ID path) or ontology synonym expansion (SNP path). Covers the full GWAS Catalog vocabulary without manual curation — the EBI ontology team maintains the synonyms |
| GWAS concurrent fetch paths + session gene cache | Primary (gene-ID) and SNP paths run via `asyncio.wait(timeout=15s)` — worst case is max(paths, 15s), not sum. Session gene cache: same gene queried for multiple traits fetches associations once and filters per trait (COX2→PTGS2/pain reuses PTGS2/inflammation fetch: 41.7s → 1.9s) |
| GWAS score saturation at 3 hits | API pagination means we may get 4 hits instead of 16 for the same target on different runs. Saturating at 3 makes the score stable to fetch boundaries |
| Biologics floor on OT score | Prevents pharmacologically-validated targets with CNV amplification or cytokine inhibition mechanisms from being underscored by OT datatypes that don't apply to their biology |
| Pydantic V2 models with `to_markdown()` | Tools output agent-readable markdown strings, never raw JSON dicts — MCP clients render output directly |
| DepMap + EFO 7-day disk cache | EFO IDs are stable across quarterly releases; DepMap CSV download (~30s) cached so warm starts are instant |
| `asyncio.Semaphore` per rate-limited API | STRING/ChEMBL: `Semaphore(2)`, Reactome: `Semaphore(3)` — prevents 429s without serializing the pipeline |
| `response_format` param on all tools | `"markdown"` default for agent consumption; `"json"` for programmatic / pipeline integration |
| `ToolSpec` registry with `use_when` | Each tool has a semantic description for dynamic selection — the workflow agent reads this to decide which tools to call for a given question |

### MCP Resource: `tool://registry`

Clients and agent frameworks can read `tool://registry` to get a structured Markdown catalogue of all 14 tools grouped by category, with `use_when` guidance for semantic retrieval. No tool call needed — useful for embedding-based tool selection or capability auditing.

---

## Development

```bash
uv sync
uv run pytest tests/ -v          # 46 unit + integration tests
uv run pytest tests/ -v --cov=genesis_bio_mcp

# Full integration test against live APIs (saves reports to examples/)
uv run python test_full.py

# Single target
uv run python test_full.py --gene BRAF --disease melanoma
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
5. Add a `ToolSpec` entry in `workflow_agent.py`'s `build_tool_registry()` with `tool_category` and `use_when`
6. Add unit tests in `tests/test_clients.py` using `respx` to mock HTTP at the transport level
7. Run `uv run ruff format . && uv run ruff check --fix .` and `uv run pytest tests/ -v`
8. Update this README

---

## Environment variables

| Variable | Default | Required for | Description |
|----------|---------|--------------|-------------|
| `ANTHROPIC_API_KEY` | — | `run_biology_workflow` | Anthropic API key for the inner Claude agent |
| `NCBI_EMAIL` | `genesis-bio-mcp@example.com` | All NCBI queries | Required for NCBI E-utils polite use |

---

## API reference

| Database | API Type | Rate Limit | Notes |
|----------|----------|------------|-------|
| UniProt | REST | Generous | `organism_id:9606 AND reviewed:true` for human Swiss-Prot |
| Open Targets | GraphQL | ~2 req/s | Requires Ensembl ID + EFO ID — resolved automatically |
| DepMap | REST (task queue) | Moderate | Celery polling → pre-signed CSV URL. 7-day disk cache. OT fallback if unreachable |
| GWAS Catalog | REST/HAL | Moderate | Concurrent gene-ID + SNP fetch paths; 15s global timeout; 24h disk cache |
| EBI OLS4 | REST | Generous | EFO ontology lookup for trait resolution. 7-day disk cache. Data queried at runtime — not redistributed |
| PubChem | REST | 5 req/s | HTTP 503 on rate limit; tenacity + Semaphore |
| ChEMBL | REST | ~1 req/s | No API key required. Two-step: target → bioactivity. `Semaphore(2)` |
| AlphaFold | REST | Generous | pLDDT ≥90 = high confidence; <70 = disordered |
| RCSB PDB | REST | Generous | Search → entry fetch per structure |
| STRING | REST | ~2 req/s | `required_score=700` for high-confidence interactions. `Semaphore(2)` |
| DGIdb | GraphQL | Generous | DGIdb 5.0; deduplicates drugs, surfaces approval status |
| ClinicalTrials.gov | REST v2 | Generous | Phase normalization: `PHASE1` → `"Phase 1"` |
| Reactome | REST | ~5 req/s | Two-step token API. `Semaphore(3)` |
| NCBI E-utils | REST | 3 req/s | Set `NCBI_EMAIL` env var |
