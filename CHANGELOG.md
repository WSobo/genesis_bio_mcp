# Changelog

All notable changes to this project will be documented in this file.

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html)

---

## [Unreleased] — v0.2.0

### Added
- `get_pathway_members` tool — enumerate all HGNC gene symbols in a named Reactome pathway by display name or stable ID; enables systematic pathway-based screening in `run_biology_workflow`
- `Settings` class (`config/settings.py`) via `pydantic-settings` — all runtime config now configurable via `GENESIS_*` environment variables or a `.env` file without code changes
- `docs/deployment.md` — full deployment guide covering env vars, Claude Desktop config, Docker setup, and production checklist
- `CONTRIBUTING.md` augmented with settings integration pattern for new clients
- **BioGRID client** — `get_biogrid_interactions` tool returning curated literature PPI records with experimental method metadata (two-hybrid, co-IP, etc.) and PubMed citations; requires `BIOGRID_ACCESS_KEY` env var (free registration at webservice.thebiogrid.org)
- `BioGRIDInteraction` and `BioGRIDInteractome` Pydantic models with `to_markdown()`
- `test_string_returns_empty_interactors_on_network_failure` — covers STRING resolve-success/network-failure edge case

### Fixed
- Reactome cache poisoning — transient network failures no longer block subsequent retries for the same gene within a session
- Reactome inline pathway parsing — switched from two-step POST→GET to using the inline `pathways` array in the `/identifiers/` response; adds trailing newline to identifier body per API contract
- `run_biology_workflow` auth error — `AsyncAnthropic()` init is now wrapped in try/except with a plain-language message directing users to set `ANTHROPIC_API_KEY` in the MCP server's env
- MCP server startup now logs a warning if `ANTHROPIC_API_KEY` is absent, surfacing the error before the first `run_biology_workflow` call

### Changed
- Server now exposes **15 tools**
- All hardcoded timeouts, cache paths, TTLs, semaphore limits, and the Claude model ID are now read from `Settings` at startup; defaults are unchanged
- Removed unused `polars` runtime dependency
- Removed stale `[project.optional-dependencies]` dev-dep block; `[dependency-groups]` is now the sole authoritative source
- Version bumped to `0.2.0`; `User-Agent` header is now dynamically read from `__version__` rather than hardcoded
- `README.md` updated with `env` block example for Claude Desktop and link to deployment guide
- `CLAUDE.md` expanded with concrete client/model/tool/test/workflow-agent patterns; dead docs/ references fixed
- `CONTRIBUTING.md` updated with current tool registration example, Step 5 (workflow agent), and updated PR checklist

---

## [0.1.0] - 2026-04-11 — Initial Alpha Release

### Summary

Genesis Bio MCP: AI-driven drug discovery via MCP tools wrapping public biomedical databases.

### Added

**Core MCP Server & Tools**

- MCP server (`server.py`) exposing 13 tools across 12 public biomedical databases
- `resolve_gene` — UniProt + NCBI alias resolution with session-level gene cache
- `get_protein_info` — UniProt Swiss-Prot curated protein metadata
- `get_target_disease_association` — Open Targets association scores
- `get_cancer_dependency` — DepMap CRISPR dependency scores with pan-essential detection
- `get_gwas_evidence` — GWAS Catalog trait matching with EFO ontology resolution
- `get_compounds` — PubChem compound search
- `get_chembl_compounds` — ChEMBL bioactivity data with pChEMBL potency scoring
- `get_protein_structure` — AlphaFold + RCSB PDB structure lookup
- `get_protein_interactome` — STRING protein interaction network
- `get_drug_history` — DGIdb + ClinicalTrials.gov drug and clinical trial history
- `get_pathway_context` — Reactome pathway enrichment
- `prioritize_target` — Composite 0–10 score across 6 evidence axes
- `compare_targets` — Multi-target ranking (2–5 targets)
- `run_biology_workflow` — Autonomous workflow agent (requires `ANTHROPIC_API_KEY`)

**Scoring Model (6 axes, summed to 10.0)**

| Axis | Max | Notes |
|---|---|---|
| Open Targets association | 3.0 | overall_score × 3 |
| DepMap CRISPR dependency | 2.0 | fraction_dependent × 2; pan-essential genes capped at 0.5 |
| GWAS evidence | 2.0 | saturates at ≥3 hits (pagination-stable) |
| Clinical / known-drug evidence | 1.5 | biologics floor: OT ≥3.25 when known_drug_score >0.9 and genetic/somatic null |
| ChEMBL potency | 1.5 | pChEMBL ≥9→1.5, ≥7→1.0, ≥5→0.5, else 0.25 |
| UniProt protein quality | 1.5 | Swiss-Prot reviewed +0.5; variant coverage up to +1.0 |

**GWAS & EFO Trait Matching**

- EFO ontology-backed trait resolution via OLS4 API (`config/efo_resolver.py`)
- Free-text disease queries resolve to canonical EFO terms (e.g. `"fat"` → obesity, `"joint inflammation"` → rheumatoid arthritis)
- 7-day disk cache for EFO term resolutions
- GWAS trait synonyms extracted to `config/trait_synonyms.py` (domain knowledge separated from HTTP client)
- `efo_uri` field added to `GwasHit` model

**Concurrency & Caching**

- `asyncio.gather` for parallel sub-queries within `prioritize_target`
- `asyncio.wait` with 15 s bound for concurrent GWAS primary + SNP fetch paths
- Session-level gene cache: eliminates redundant cross-trait fetches (COX2/pain: 41.7 s → 1.9 s)
- 7-day disk cache for DepMap CRISPR data
- 24-hour disk cache for GWAS associations (timeout resilience)
- Single shared `httpx.AsyncClient` for connection pooling
- `asyncio.Semaphore` per rate-limited API (prevents 429s)

**Output & API Design**

- Pydantic V2 models with `to_markdown()` for all structured outputs
- `response_format` param (`markdown` default | `json` for programmatic use) on all tools
- All MCP tools output strictly formatted Markdown strings (never raw dicts)
- `safe_call` wrapper on all coroutines — single API failure does not crash the pipeline
- All tools annotated: `readOnlyHint`, `destructiveHint`, `idempotentHint`, `openWorldHint`
- Tool naming convention: `{service}_{action}_{resource}` snake_case

**CI/CD & Developer Experience**

- GitHub Actions CI pipeline (lint + test on every push and PR)
- 60% line coverage threshold enforced via `pytest-cov`
- `uv`-based package management (no pip/conda)
- `ruff` for formatting and linting
- `ToolSpec` registry with `use_when` descriptions for semantic tool selection by the workflow agent

**Example Reports**

14 gene–disease example reports (JSON + Markdown) covering:

| Gene | Indication | Score | Tier |
|---|---|---|---|
| PCSK9 | hypercholesterolemia | 9.1 | High |
| PTGS2 / COX2 | pain | 7.8 | High |
| TNF | rheumatoid arthritis | 7.8 | High |
| PTGS2 | inflammation | 7.7 | High |
| HER2 / ERBB2 | breast cancer | 7.6 | High |
| EGFR | non-small cell lung carcinoma | 7.5 | High |
| BRAF | melanoma | 7.2 | High |
| KRAS | pancreatic cancer | 4.9 | Medium |
| FTO | obesity | 4.3 | Medium |
| CD274 | melanoma | 4.1 | Medium |
| TP53 | squamous cell carcinoma | 3.6 | Low |
| TP53 / p53 | lung cancer | 3.4 | Low |

**Tests**

- 46 unit + integration tests across 3 test modules (`test_clients.py`, `test_e2e.py`, `test_workflow_agent.py`)
- Mocked HTTP via `respx` at transport level
- Coverage: gene resolution, all tools, edge cases, EFO trait matching, GWAS cache correctness, session cache poisoning prevention

### Development Status

Alpha — API subject to change.

---

[0.1.0]: https://github.com/WSobo/genesis-bio-mcp/releases/tag/v0.1.0
