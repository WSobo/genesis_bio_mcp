# Changelog

All notable changes to this project will be documented in this file.

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html)

---

## [0.3.4] — 2026-04-25

Patch release. Fixes five bugs caught by a fourth-round live MCP smoke test
(EGFR, GLP1R, CALCA, TP53, ABL1, MT-ND1, RAS family, NSCLC, PDAC). The most
impactful is **Bug M.1**: the `prioritize_target` Markdown report was
recomputing every score axis inline using stale formulas, producing a
breakdown table whose numbers drifted from the canonical `_compute_score`
output. That display recompute was also masking the v0.3.3 fixes (functional
vs binding split, fallback suppression, monogenic credit) — the actual
priority score reflected them, but the table didn't.

### Fixed

- **Display drift in `prioritize_target` scoring breakdown (Bug M.1).**
  `TargetPrioritizationReport.to_markdown` was a duplicate of
  `_compute_score` with stale logic — GWAS capped at 10 instead of 3, no
  monogenic credit, no fallback suppression, no functional/binding split.
  Now renders directly from `score_breakdown` (the canonical
  `ScoreBreakdown` produced by `_compute_score`), with descriptive notes
  derived from the source models. Eliminates two parallel scoring code
  paths and surfaces the v0.3.2/v0.3.3 fixes in the displayed table for
  the first time. Also adds a "Clinical / known-drug" row (was previously
  absorbed into OT) and explanatory notes for fallback / off-trait
  suppression so a 0.00 GWAS score isn't mistaken for "no data."
- **GWAS axis credited off-trait direct-match returns (Bug M refined).**
  v0.3.3 zeroed GWAS credit when the fallback path explicitly marked hits
  as off-trait via the `trait_query` sentinel. But TP53/Li-Fraumeni
  produced direct-match GWAS returns whose trait labels were
  sex-hormone-binding-globulin and height — completely unrelated to
  Li-Fraumeni, but `_process_for_trait` let them through. Added a
  trait-relevance gate in `_compute_score`: token-set Jaccard between the
  indication and each hit's trait label, with a `_GWAS_TRAIT_RELEVANCE_MIN`
  threshold (0.15) below which the GWAS axis is zeroed. Catches
  direct-match returns where the queried indication substring-matched
  something incidental in GWAS Catalog. Calibrated so PCSK9/CHD (~0.50)
  and EGFR/NSCLC (~0.20) keep credit; TP53/sex-hormone (~0.0) gets zeroed.
- **OT silent-data-loss on bare acronyms (Bug P).** `EGFR + "NSCLC"`
  returned 0 OT evidence because OT's autocomplete-style `search` doesn't
  index abbreviations. Different failure mode from Bug J (parens
  contamination) — Bug P is the bare acronym alone. Added
  `_ACRONYM_EXPANSIONS` table with ~50 common indication abbreviations
  (NSCLC, PDAC, T2DM, AD, PD, CFTR-relevant: CF, AML, CML, COPD, NASH/MASH,
  HFpEF/HFrEF, etc.) keyed off all-caps tokens in the input. The
  expansion is added as a fallback variant to `_normalize_indication_variants`
  so the literal acronym is still tried first (cache-friendly).
- **Sub-nM IC50 rendered as "0.0 nM" (Bug Q).** Calcitonin gene-related
  peptide / GLP1R agonists hit pChEMBL ≈ 11 (sub-pM territory) which the
  format string rounded to "0.0 nM" — looked like a missing-data error
  rather than ultrahigh potency. New `_format_ic50_nm()` helper renders
  sub-nM in pM, 1–1000 nM in nM, beyond in µM. Used by both the
  ChEMBLCompounds table and the prose summary in `target_prioritization`.
- **Reactome cross-species pathway leakage (Bug R).** The AnalysisService
  doesn't have a species filter on the request; for human gene queries it
  could return canine (`R-CFA-*`), mouse (`R-MMU-*`), or rat (`R-RNO-*`)
  pathways when the gene mapped across species silos. GLP1R surfaced an
  `R-CFA-381676` row in its top pathway list. `_parse_pathways` now
  drops anything whose stable ID doesn't start with `R-HSA-` (also drops
  empty-stId entries which weren't actionable anyway).

### Verified — no fix needed

- **Bug S** — `run_biology_workflow` already handles missing
  `ANTHROPIC_API_KEY` gracefully: server logs a warning at startup
  (`server.py:110`) and the agent itself returns a clear
  "Set ANTHROPIC_API_KEY in claude_desktop_config.json" message rather
  than a stack trace. The user explicitly OK'd this behavior.

### Added

- `_format_ic50_nm()` helper in `models.py` for unit-aware IC50 rendering.
- `_ACRONYM_EXPANSIONS` table in `clients/open_targets.py` (~50 entries
  spanning oncology, cardiometabolic, neuroscience, immunology).
- `_GWAS_TRAIT_RELEVANCE_MIN` constant + `_max_trait_relevance` /
  `_tokenize_for_relevance` helpers in `tools/target_prioritization.py`
  with documented tuning rationale.
- 5 regression tests in `tests/test_clients.py` (one per bug, including a
  display-vs-breakdown drift test that asserts the rendered table matches
  the canonical `ScoreBreakdown` numbers).
- **214/214 tests pass.**

### Still deferred to v0.4.0

- **Bug N** — `compare_targets` DepMap column shows blank for pan-essential
  rows (cosmetic).
- **Bug O** — Compounds column conflates PubChem and ChEMBL counts (cosmetic).
- **Bug H** — GWAS empty-result latency optimization.
- **EFO → UBERON** ontology-backed indication-to-tissue mapping.
- **5-target live-API regression check** + **CDR developability** + **Foldseek**
  + **per-residue pLDDT** + **dedicated MaveDB tool**.

---

## [0.3.3] — 2026-04-25

Patch release. Fixes six bugs caught by a third-round live MCP smoke test
(EGFR, ABL1, GLP1R, MT-ND1, BCR-ABL1, TP53/Li-Fraumeni, KRAS/NRAS/HRAS/MRAS
in PDAC). Includes one **silent-data-loss critical** in EFO matching and
two scoring-rubric corrections that were affecting target rankings on
real-world queries. Also folds in the chemical-matter assay-type weighting
that was deferred to v0.4.0.

### Fixed

- **OpenFDA off-target leak (Bug F.1).** v0.3.2 restricted FAERS enrichment
  to direct-engagement DGIdb interaction types but left an "untyped
  fallback" for newer approvals DGIdb hadn't categorized yet. That fallback
  let in spurious co-mentions: ACYCLOVIR and ATP appeared in ABL1's safety
  panel, and ATEZOLIZUMAB / BEVACIZUMAB still leaked into EGFR's. Removed
  the untyped fallback entirely — DGIdb leaves `interaction_type=None` for
  exactly the spurious associations we don't want, so accepting untyped
  defeats the filter. Trade-off: very-recent approvals not yet typed by
  DGIdb won't get FAERS coverage; this is the correct precision/recall
  trade-off for a safety panel.
- **OpenFDA selection by FAERS volume, not alphabetical (Bug F.2).** v0.3.2
  sorted candidates by `(phase desc, name asc)`, so for ABL1 (six phase-4
  inhibitors) ASCIMINIB beat IMATINIB on alphabetical sort despite IMATINIB
  having decades more FAERS data. Now queries a wider candidate pool (10
  vs 5), pre-ranks by `(phase desc, source-count desc, name asc)`, then
  re-ranks the returned signals by `total_reports` and keeps the top 5
  for display. Result: established drugs surface ahead of recent approvals.
- **EFO matcher fuzzy-substring silent-data-loss (Bug K).** Same class as
  Bug J. `"type-2 diabetes mellitus (T2DM)"` returned `MONDO_0010802`
  (pancreatic-hypoplasia-diabetes-congenital-heart-disease syndrome) — a
  rare congenital syndrome that fuzzy-matched on the word "diabetes". GLP1R
  (the most validated T2D target) scored 0.767 for the wrong disease.
  Resolver now scores every variant's hits by token-set Jaccard against
  *both* the variant and the original input, picks the highest-scoring
  match across all variants, and rejects matches below
  `_DISEASE_MATCH_THRESHOLD` (0.5). The pancreatic-hypoplasia hit scores
  ~0.13 and is rejected; the real T2D entry scores 1.0 and wins.
- **Pan-essential cap missed mitochondrial-encoded genes (Bug L).** MT-ND1
  scored 2.0/2.0 cancer-dependency credit despite being a mitochondrial
  complex I subunit (you can't drug it without destroying every
  mitochondrion). DepMap's `common_essential` flag doesn't fire for genes
  outside the routine CRISPR screen. Added a defensive heuristic: if
  ≥95% of measured cell lines are dependent, treat as pan-essential
  regardless of the `common_essential` flag. Applies to both the cache
  path and the OT-proxy fallback path.
- **GWAS fallback over-credit (Bug M).** v0.3.2's Bug D fix surfaced
  unfiltered top gene-level GWAS hits when no exact-trait match was found.
  But `_compute_score` then awarded full GWAS credit for those hits even
  though they were under unrelated traits. Result: TP53/Li-Fraumeni earned
  GWAS credit from sex-hormone-binding-globulin and height GWAS hits at
  the TP53 locus; MRAS outranked HRAS for PDAC because off-trait hits
  inflated its score. Fix: detect the `"no exact-trait match"` sentinel
  in `gwas_ev.trait_query` and zero the GWAS axis. Hits still appear in
  the report for context — they just don't score.
- **Chemical-matter scoring didn't weight assay type (Bug I, v0.4.0
  pull-forward).** MYC scored 1.5/1.5 from binding-only assays (MYC-MAX
  heterodimerization, G-quadruplex binders) despite no functional /
  cellular activity — those binders virtually never translate. Added
  `best_pchembl_functional` and `best_pchembl_binding` fields to
  `ChEMBLCompounds`, populated from `assay_type` (F vs B) in the ChEMBL
  client. `_compute_score` now prefers functional potency at full credit
  (1.5 / 1.0 / 0.5 / 0.25 by tier) and discounts binding-only by ~33%
  (1.0 / 0.7 / 0.35 / 0.2). Targets with strong functional readouts
  (kinases) keep their scores; binding-only undruggables (MYC,
  β-catenin) lose ~0.5 credit on the chemical-matter axis.

### Added

- `_DISEASE_MATCH_THRESHOLD = 0.5` and `_name_match_score()` helper in
  `clients/open_targets.py` with documented scoring (substring rule for
  acronyms, Jaccard otherwise).
- `_SAFETY_LOOKUP_POOL = 10` in `tools/target_prioritization.py` —
  candidate pool for the new two-pass FAERS ranking.
- `best_pchembl_functional` and `best_pchembl_binding` fields on
  `ChEMBLCompounds` model with documented tuning rationale.
- 6 regression tests in `tests/test_clients.py` (one per bug).
- 2 existing tests updated to match intentional behavior changes
  (DepMap pan-essential cap firing in mocked 100% dependency case;
  attach_safety_signals strict-filter dropping untyped drugs).
- **209/209 tests pass.**

### Still deferred to v0.4.0

- **Bug N** — `compare_targets` DepMap column shows blank for pan-essential
  rows; cosmetic, no scoring impact.
- **Bug O** — Compounds column conflates PubChem and ChEMBL counts;
  cosmetic clarity issue in `compare_targets` table.
- **Bug H** — GWAS empty-result latency optimization.
- **EFO → UBERON** ontology-backed indication-to-tissue mapping
  (replaces hardcoded `config/indication_tissue_map.py`).
- **5-target live-API regression check** (BRAF/melanoma, EGFR/NSCLC,
  PCSK9/hypercholesterolemia, TNF/RA, KRAS/pancreatic).
- **CDR developability** on `get_antibody_structures`, **Foldseek**
  structural homologs, **per-residue pLDDT bands**, **dedicated
  per-variant MaveDB tool**.

---

## [0.3.2] — 2026-04-25

Patch release. Fixes five bugs caught by a second-round live MCP smoke test
(PCSK9, CFTR, MYC, EGFR, EVOLOCUMAB, EGFR/NSCLC, INS, MIR21, OCT4, MLL),
including one **silent-data-loss critical** in Open Targets and one v0.3.1
fix that turned out to be insufficient (GTEx).

### Fixed

- **GTEx still returned empty TPM payloads (Bug B carryover from v0.3.1).**
  v0.3.1 fixed the unversioned-ID case by routing through GTEx's
  `/reference/gene` endpoint. But that endpoint defaults to
  `gencodeVersion=v26`, which returns IDs the current expression dataset
  (`gtex_v10`, keyed on **GENCODE v39**) can't find — the lookup succeeded
  but the expression endpoint silently returned `data: []` again. Pin
  `gencodeVersion=v39` so the IDs match what `gtex_v10` indexes against.
  Verified end-to-end against PCSK9 (liver TPM ~5–10) and INS (pancreas
  TPM ~10000+).
- **Open Targets indication parsing — silent data loss on parentheses.**
  `EGFR + "non-small-cell lung cancer (NSCLC)"` returned 0 OT evidence and
  scored 3.4 LOW PRIORITY for what is the textbook precision-oncology
  success story. Root cause: OT's `search` GraphQL is autocomplete-style
  and brittle — both `"non-small-cell lung cancer (NSCLC)"` AND
  `"non-small-cell lung cancer"` return zero hits, but `"NSCLC"` alone
  resolves to `EFO_0003060`. Added `_normalize_indication_variants()` that
  tries the literal first (cache-friendly), then the parenthetical-stripped
  form, then the bare abbreviation, then the hyphen-normalized form,
  logging which variant won. Worst class of bug closed: confident-looking
  wrong answer on a real target.
- **GWAS trait-name mismatch (Bug D carryover).** When the queried
  indication doesn't match any GWAS Catalog trait label (PCSK9 has dozens
  of LDL/lipid associations, none labelled `"familial hypercholesterolemia"`),
  the trait filter returned None and the GWAS axis zeroed out. Added a
  fallback in `GwasClient.get_evidence` that surfaces the strongest
  unfiltered gene-level associations with a sentinel `trait_query` so the
  score reflects "gene IS GWAS-implicated, just not under this exact label."
- **Composite scoring penalized monogenic diseases.** `CFTR + cystic
  fibrosis` scored 6.7 MEDIUM despite being the most-validated CF target —
  the GWAS axis was 0 because GWAS Catalog doesn't study Mendelian
  diseases. Added a monogenic credit: when OT's
  `genetic_association_score > OT_GENETIC_MONOGENIC_THRESHOLD` (0.7) and
  GWAS returns nothing, award the full 2.0 GWAS axis. Polygenic targets
  (FTO, PCSK9) sit below the threshold and are unaffected.
- **FAERS pulled non-target drugs from DGIdb co-mentions.** `get_drug_history`
  for EGFR returned safety panels including atezolizumab (anti-PD-L1) and
  bevacizumab (anti-VEGF) because DGIdb associates them with EGFR via
  trial co-administration. Restricted `attach_safety_signals()` to drugs
  whose `interaction_type` is in `_DIRECT_TYPES`
  (inhibitor / antagonist / blocker / agonist / modulator / binder), with
  untyped approved drugs retained as fallback for newer approvals DGIdb
  hasn't categorized yet.

### Added

- 6 regression tests in `tests/test_clients.py` — one per bug (B, J, D, F,
  G) plus a unit test for `_normalize_indication_variants`. **203/203 tests
  pass.**
- `OT_GENETIC_MONOGENIC_THRESHOLD = 0.7` constant in
  `tools/target_prioritization.py` with documented rationale and tuning
  range, alongside the existing `OT_CLINICALLY_VALIDATED_FLOOR`.

### Deferred to v0.4.0

- **Bug H** — GWAS empty-result latency (~10s on no-hit cases). Not a real
  bug; the SNP-fetch + cascading association fetches have a hard 15s ceiling.
  Bug D's unfiltered fallback partially masks this since it returns earlier.
- **Bug I** — chemical-matter score doesn't weight assay type. MYC scores
  1.5/1.5 from binding-only assays despite being undruggable. Real
  correctness gap, but fixing it requires a new `best_pchembl_assay_type`
  field on `ChEMBLCompounds` and recalibration against the 14-target
  benchmark to confirm no tier flips. v0.4.0 scoring overhaul, not a patch.

---

## [0.3.1] — 2026-04-25

Patch release. Fixes four bugs caught by the post-v0.3.0 smoke run against
PCSK9, BRAF, EGFR, and HMGCR. All four were present in v0.3.0; none change
the public tool surface.

### Fixed

- **GTEx returned empty TPM payload for every gene** — the expression
  endpoint silently returns `data: []` when given an unversioned GENCODE
  ID (`ENSG…`). v0.3.0 resolved through Ensembl which returns the bare ID;
  the version suffix (e.g. `.11`) is GTEx-specific and mandatory. Resolution
  now goes through GTEx's own `/api/v2/reference/gene` endpoint, which
  returns the GENCODE version GTEx itself indexes against. EnsemblClient is
  kept as a fallback for symbols GTEx doesn't know.
- **OpenFDA safety section silently dropped for biologics + salt forms** —
  two compounding bugs:
  1. The Lucene query joined three `field:value` clauses with `+`, which
     Lucene parses as **MUST** (AND), so a record had to match the drug
     name in `medicinalproduct`, `openfda.generic_name`, *and*
     `openfda.brand_name` simultaneously. Worked for clean INNs like
     `fluvastatin`; failed for `alirocumab`, `evolocumab`, and TKIs whose
     spelling varies across fields. Replaced `+` with explicit `OR` and
     added `openfda.substance_name` as a fourth match field.
  2. DGIdb returns drug names with salt and pharmaceutical-form suffixes
     (`ATORVASTATIN CALCIUM TRIHYDRATE`, `DACOMITINIB ANHYDROUS`). OpenFDA's
     `openfda.*` fields don't index these qualifiers. Added
     `_normalize_drug_name()` which iteratively strips a list of common
     salts, hydrates, and counterions before the lookup.
- **UniProt resolver picked the wrong gene for some symbols** — `gene_exact:ALB`
  returns `FBF1` ranked first in UniProt's relevance ordering; the resolver
  blindly took `results[0]`. Fetch widened from 1 to 5 results, with a new
  `_pick_exact_gene_match()` helper that prefers entries whose primary
  `geneName` matches the query before falling back to first hit. Affected
  any single-gene tool through `_resolve_symbol`.
- **ChEMBL `assay_organism` always rendered empty** — the activity row
  exposes the species via `target_organism` (populated); `assay_organism`
  is almost always null in the live API. Parser now reads `target_organism`
  first, falling back to `assay_organism` only when the assay system differs
  from the target species (e.g. human protein expressed in insect cells).

### Added

- 5 regression tests in `tests/test_clients.py` — one per bug, plus a unit
  test for `_normalize_drug_name`. **197/197 tests pass.**

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
