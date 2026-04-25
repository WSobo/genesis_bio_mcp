# Changelog

All notable changes to this project will be documented in this file.

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html)

---

## [0.3.0] — 2026-04-22

Major feature release. Closes the genomics coordinate gap, adds tissue/protein expression evidence, introduces a post-market drug-safety layer, and deepens ChEMBL assay context. Two new MCP tools (`get_variant_consequences`, `get_tissue_expression`, `get_protein_atlas`) plus meaningful depth fixes to `get_variant_effects`, `get_drug_history`, `get_chembl_compounds`, and `prioritize_target`. No breaking changes to tool names or argument shapes — all additions are additive on the response side.

### Added

- **Ensembl REST client + VEP (M1)** — `clients/ensembl.py` wraps `/lookup/symbol`, `/overlap/region`, `/vep/human/hgvs` with a 5-req/s semaphore. `EnsemblGene`, `VEPConsequence`, `TranscriptInfo` models with `to_markdown`. New MCP tool `get_variant_consequences` returns splice/UTR/regulatory overlap and VEP's own SIFT/PolyPhen — complementary to the dbNSFP values from MyVariant. Canonical-transcript by default; `include_all_transcripts` opt-in for full coverage.
- **GTEx + Human Protein Atlas clients (M2)** — `clients/gtex.py` fetches median TPM per tissue (resolves HGNC → GENCODE via the shared EnsemblClient). `clients/hpa.py` parses the HPA bulk-download JSON for tissue-specificity category, subcellular localization, and pathology prognostics. New MCP tools `get_tissue_expression` and `get_protein_atlas`. Both share the 7-day disk-cache pattern from SAbDab. New `config/indication_tissue_map.py` with a hardcoded top-20 therapeutic-area → tissue mapping (interim; v0.3.1 plans EFO → UBERON resolution).
- **OpenFDA drug-safety client (M3)** — `clients/openfda.py` wraps FAERS `/drug/event.json`, structured label `/drug/label.json` (boxed warnings), and `/drug/enforcement.json` (recalls). Fans out all three sub-queries in parallel with a 2-concurrency semaphore. 7-day disk cache. Optional `OPENFDA_API_KEY` env var lifts the 240 req/min / 1000 req/day free quota. `AdverseEventCount`, `DrugRecall`, `DrugSafetySignal` models, all carrying a permanent disclaimer field so FAERS counts are never rendered without the regulatory caveat. No new MCP tool — `get_drug_history` and extended-mode `prioritize_target` now populate `.safety` on the top-5 approved drugs by phase via a new `tools/target_prioritization.attach_safety_signals()` helper that isolates per-drug failures.
- **Expression axis in priority scoring (M5)** — `ScoreBreakdown` gains an `expression` field (max 1.0). `_compute_score` accepts an optional `protein_atlas` argument and applies the HPA-derived tissue-specificity bonus via a new `_EXPRESSION_BY_CATEGORY` table (`Tissue enriched` = 1.0 → `Not detected` = 0.0). `prioritize_target` fetches the HPA report in extended mode and passes it through; the 10.0 score cap absorbs the new axis so existing tiers do not regress.
- Test coverage added this release: 4 Ensembl + VEP tests, 3 GTEx tests, 3 HPA tests, 4 OpenFDA tests, 4 ChEMBL tests (previously uncovered), 9 scoring-axis parametrized tests. **184/184 tests pass on `feat/v0.3.0`.**

### Changed

- **ChEMBL assay context depth (M4)** — `ChEMBLActivity` gains `assay_type` (B/F/A/T/P), `assay_organism`, `assay_cell_type`, `bao_format`, and `confidence_score`. All five fields are already in the ChEMBL activity response; they were dropped before and are now parsed. `ChEMBLCompounds.to_markdown` surfaces an **Assay mix** summary (e.g. `12 binding, 5 functional (3 cell-based)`), flags non-human organisms, and flags `confidence_score < 9`. The per-row table gains Assay and Organism columns; functional rows render cell type inline (e.g. `F (A375)`). Biochemical and cell-based potency numbers are no longer visually identical.
- **`get_variant_effects`** — third parallel task added: VEP consequences via `EnsemblClient.get_vep_consequences()`. `VariantEffects` gains a `vep_consequences` field; the aggregator no longer depends on a gnomAD hit to return a non-empty result.
- **`prioritize_target` signature** — new keyword `hpa: HPAClient | None = None` (parallels `reactome`, `openfda`). Extended-mode fetches now include HPA so the expression axis can score it. Fan-out happens before `_compute_score` so the HPA report feeds the scoring table.
- Version bumped to `0.3.0`; `__init__.py` drift from `0.2.1` also fixed.

### Fixed

- **MyVariant.info vs. gnomAD routing** — AlphaMissense removed from the gnomAD GraphQL path (the v4 schema does not expose it); AlphaMissense, ClinVar, REVEL, CADD, SIFT, PolyPhen, and gnomAD AF are now sourced exclusively from MyVariant/dbNSFP.
- **IEDB transport** — async-ticket pattern via the HTTPS NextGen Tools host (`api-nextgen-tools.iedb.org`); the deprecated `tools-cluster-interface.iedb.org` direct POST path is removed.

### Deferred to v0.3.1 / v0.4.0

- 5-target live-API regression check (BRAF/melanoma, EGFR/NSCLC, PCSK9/hypercholesterolemia, TNF/RA, KRAS/pancreatic) — scoring invariants covered at unit level in M5.
- EFO → UBERON ontology-backed indication-to-tissue mapping (hardcoded top-20 table lives in `config/indication_tissue_map.py` as an interim solution).
- CDR developability on `get_antibody_structures`, Foldseek-based structural homologs, per-residue pLDDT bands, dedicated per-variant MaveDB tool.
- ProteomicsDB / CPTAC, OMIM / ClinGen, patent-landscape, AlphaFold Multimer, ESM Metagenomic Atlas.

---

## [0.2.4] — 2026-04-17

Polish release addressing six issues surfaced by a JAK2 end-to-end evaluation. All fixes use dynamic, structural solutions (EFO URIs, activity_outcome fallback, token-prefix dedup) rather than hardcoded per-trait/per-drug vocabulary.

### Fixed

- **GWAS trait matching** — `EFOResolver.resolve()` now hierarchy-expands each resolved term via OLS4's `allChildrenOf` + `ancestorsOf` filters and stores the URI set on `EFOTerm.related_uris`. `filter_by_trait` matches hits against the expanded set, so `"polycythemia vera"` catches JAK2 studies tagged with `"myeloproliferative neoplasm"` (direct EFO parent), and `"myeloproliferative"` catches studies tagged with specific subtypes (EFO descendants). No new hardcoded synonyms; the fix generalizes to any indication EFO covers.
- **`get_compounds` blank Activity column** — `Compounds.to_markdown` now falls back to `activity_outcome` ("Active") when PubChem's concise endpoint omits the Activity Name cell for some assay rows.
- **`get_antibody_structures` "0.0 nM" summary** — `_parse_float` now rejects zero, negative, and non-finite values; the summary filter requires `> 0`, so the "Best measured affinity" insight is omitted rather than rendered as `0.0 nM` when affinities are unreported.
- **`get_drug_history` salt-form duplicates** — new `_collapse_salt_forms()` merges DGIdb records whose `drug_name`'s first whitespace token matches a shorter single-token record's full name (pharma salt convention, e.g. `FILGOTINIB` + `FILGOTINIB MALEATE` → one row). No hardcoded salt vocabulary.
- **`get_pathway_context` duplicate pathway names** — `_parse_pathways` now runs a second-pass dedup on `display_name` (case-insensitive), keeping the row with the smallest p-value. Reactome stable IDs are appended to the rendered pathway name so any residual duplicates are visually distinguishable.

### Added

- **`compare_targets` per-row score breakdown** — new `ScoreBreakdown` model capturing contributions from each of the six scoring axes (OT, DepMap, GWAS, known-drug, chem matter, protein). `_compute_score` returns the breakdown; it's stored on `TargetPrioritizationReport` and `TargetComparisonRow`, and `ComparisonReport.to_markdown` renders a compact per-row line (`OT 2.3 · Dep 1.2 · GWAS 0.0 · Drug 1.1 · Chem 1.5 · Prot 0.6`). Makes rankings auditable — a target with the highest OT score that ranks below peers now shows exactly which axes its peers won on.
- Eight targeted regression tests in `tests/test_clients.py` covering each fix above (EFO hierarchy expansion, URI-set matching, salt-form merging, pathway name dedup, compounds fallback, antibody affinity guard, score breakdown invariant).

### Changed

- Version bumped to `0.2.4`; `EFOTerm` gains a `related_uris: list[str]` field (default empty, persisted in the 7-day EFO disk cache).

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
