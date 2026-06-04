# Roadmap

Current release: **v0.7.0** — 42 tools across 24 public biomedical data sources (plus local
RDKit cheminformatics, the UMA-Inverse model service, and an optional embedding-backed corpus),
with `prioritize_target` / `compare_targets` orchestration (now a transparent per-axis **evidence
profile** rather than a composite score), the `run_biology_workflow` agent, the agent-contract
provenance/typed-error envelope on all JSON output, and the v0.6.0 hybrid-retrieval `corpus_*` tools
(Postgres + pgvector).

## Evidence-driven direction

Roadmap priorities are now set by the **evaluation harness** ([docs/eval.md](eval.md)), not by
hand-poking. Its baseline run answers 23/24 tasks across all eight tool categories and flags four
**capability gaps** no current tool can address — **ADMET/PK prediction**, **kinome-wide selectivity**,
**patents / freedom-to-operate**, and **retrosynthesis**. These are the candidate next builds; the gap
report is the artifact that decides which leads.

## North star

> Evolve genesis-bio-mcp from a target-evidence **retrieval** server into an AI-native,
> production-grade discovery **platform** — adding small-molecule/cheminformatics depth,
> *computed* methods (not just database lookups), and the observability, provenance, and
> agent-contract guarantees that both human scientists and autonomous agents need to trust it.

This direction is shaped by what AI-native drug-discovery orgs (e.g. Eli Lilly's Data Foundry /
Small Molecule Discovery, Genesis Molecular AI) actually ask of their tooling — which is
**not "more public databases"** but reliable, governed, agent-consumable *methods*:

| Industry pillar | Implication for this MCP |
|---|---|
| **Methods4Insight** — analytical/computational methods | Tools that *compute* on structures, not only retrieve |
| **Architecture4Insight** — data infra + scientific software | Clean APIs; a seam for *proprietary* data behind the same interface |
| **Automation & Scale4Insight** — agentic workflows | Batch ops, deterministic contracts; "serve humans *and* autonomous agents" |
| **Preparedness4Insight** — governance/readiness | Provenance, versioning, reproducibility, audit trail |
| MLOps practice | Model serving, observability guardrails, response-time guarantees, typed error handling |

The biggest current gap: the server is ~90% target/biology evidence with only thin
small-molecule lookups, while the highest-value industrial use case (small-molecule
discovery) needs cheminformatics. **v0.5.0 closes that gap.**

---

## v0.5.0 — Cheminformatics core + a served inverse-folding model + agent/ops hardening

**Theme:** small-molecule capability (RDKit, offline & deterministic), a first **served ML model**
(UMA-Inverse, structure → sequence), and the provenance/observability/agent-contract slice that
makes all 31+ tools enterprise- and agent-grade. Together these cover all three flavors of tool a
modern discovery stack needs: **retrieve** (databases), **compute** (RDKit), and **infer** (a
deployed model).

**Dependency budget (decided):** **RDKit only** for the local cheminformatics. RDKit (BSD,
pure-Python wheel) unlocks offline, deterministic, rate-limit-free chemistry. The UMA-Inverse
integration (M4) adds **no heavy dependency** — it's a thin REST client to an already-deployed
service. Model-based ADMET and docking (which *would* add heavy ML deps) remain deferred to v0.6.

### Track 1 — Cheminformatics & computed methods *(flagship; Methods4Insight)*

| Milestone | Tool | Description |
|---|---|---|
| **M1** | `compute_molecular_properties` | SMILES → MW, logP, TPSA, HBD/HBA, rotatable bonds, aromatic rings, **QED**, **Lipinski/Veber/PAINS flags**, **Bemis-Murcko scaffold**. Pure RDKit, fast, deterministic. |
| **M2** | `standardize_structure` | Salt stripping, tautomer/charge normalization, canonical SMILES + **InChIKey**. Data-readiness primitive other tools build on (`Preparedness4Insight`). |
| **M3** | `search_similar_compounds` | Morgan-fingerprint **Tanimoto similarity** + **substructure (SMARTS)** search over ChEMBL/PubChem actives for a target or query molecule. |
| **(stretch)** | `analyze_sar` | Matched-molecular-pairs / activity-cliff summary over a ChEMBL target's compound set. |

This turns the server into a **closed loop**: target validation → known chemical matter →
property/scaffold analysis → similar-compound expansion.

### Track 2 — Inverse folding & model serving *(M4; Methods4Insight + serving/deployment)*

Integrate **[UMA-Inverse](https://huggingface.co/spaces/WSobo/uma-inverse)** — a deployed,
containerized, monitored inverse-folding service (structure → sequence) — as first-class
genesis-bio-mcp tools. This adds a *served ML model* to the stack (not just a database or a
local compute), and is the clearest embodiment of the industry "model serving + observability +
API layer" pillar: UMA-Inverse already ships FastAPI, Prometheus `/metrics`, `structlog` JSON
logs, `/health`, OpenAPI `/docs`, a Gradio UI, and a Docker image (HF Spaces CPU Basic, ~50 ms
to design a 46-residue protein).

| Milestone | Tool | Description |
|---|---|---|
| **M4** | `design_sequence_for_structure` | `POST /design` — redesign a backbone's sequence. Input: full PDB text (+ optional `ligand`, `temperature`). Returns designed `sequences`, per-residue + mean confidence, residue count, `inference_ms`, `request_id`. |
| **M4** | `score_structure` | `POST /score` — score a sequence against a structure (native by default). Returns per-residue `log_prob`/`prob`, the model's preferred residue (`top_aa`/`top_prob`), overall `perplexity`, and `recovery` — surfaced as a **candidate-mutation table** (residues the model would change). |

**Why it closes the loop:** an agent retrieves a structure (`get_protein_structure` /
`get_structure_confidence` / AlphaFold), **scores** it to flag suboptimal residues, then
**redesigns** it — end-to-end protein engineering against a live model, complementing the
existing `get_structure_confidence` / `get_cdr_developability` / `get_dms_*` tools.

**Implementation notes:** a new `UMAInverseClient` wrapping the REST service (base URL from
`UMA_API_URL`, default the live Space `https://wsobo-uma-inverse.hf.space`); standard never-raise
handling for the documented statuses (`413` over the residue cap, `422` bad body, `504` timeout)
with a graceful message if the endpoint is unreachable; `SequenceDesign` / `StructureScore`
models with `to_markdown`; respx-mocked tests for both endpoints. Tools annotated
`openWorldHint=True` (calls an external model service). The live Space caps inputs at 256
residues (`UMA_MAX_RESIDUES`); larger workloads target a GPU endpoint. See the service's
`docs/DEPLOY.md` / `docs/SERVING_NOTES.md`.

### Track 3 — Agent contract & provenance *(Preparedness4Insight + agent-readiness)*

Touches all tools, so do it once, early.

- **M5 — Provenance envelope.** Every `response_format="json"` result carries a consistent
  block: `source`, **upstream DB/dataset version**, `retrieved_at`, `confidence`, and the
  resolved query. (Today version/release is surfaced ad hoc — DepMap release, OT v4,
  AlphaFold model version; formalize into one schema.)
- **M6 — Typed error taxonomy.** Replace prose error strings with machine-readable statuses
  (`NotFound` / `InvalidInput` / `RateLimited` / `UpstreamUnavailable`) so agents branch
  deterministically instead of string-matching.

### Track 4 — Observability & response-time guarantees *(MLOps practice)*

- **M7 — `health://status` resource** reporting per-upstream reachability + cache state, and a
  metrics surface (per-tool call count / p50–p95 latency / error rate). Builds on the existing
  per-call `api_latency_s`.
- **M8 — Batch tools** for agent fan-out: `batch_resolve_genes`, `batch_compute_properties`
  (and `batch_prioritize_targets` if cost allows) — N round-trips → 1.

### Suggested sequence
1. **M1** + **M2** (RDKit foundation — unblocks all chem) — *M1 shipped*
2. **M3** (similarity/substructure — most-requested chem capability)
3. **M4** (UMA-Inverse integration — closes the structure→score→redesign loop; the marquee
   "served model" capability)
4. **M5** + **M6** (agent contract; cross-cutting, do once) — *shipped: M5 provenance envelope, M6 typed error taxonomy*
5. **M7** + **M8** (ops/agent-readiness) — *shipped: M7 `health://status`, M8 batch tools*
6. `analyze_sar` stretch if time permits

Each milestone follows the established recipe (model + client/tool + `run_biology_workflow`
registration + tests + count-sync + CHANGELOG) and ships as its own atomic PR.

---

## v0.6.0 — Hybrid retrieval (embedding-backed corpus) *(committed)*

**Theme:** add a second retrieval pattern alongside the live-API tools — a persistent,
indexed **PostgreSQL + pgvector** store over a *curated bioactivity corpus*, exposed through
new `corpus_*` MCP tools that do **hybrid retrieval** (relational SQL filters + similarity
search) and return results with persistent IDs, source/version provenance, and honestly-
separated confidence signals. The existing 38 tools proxy live APIs and stay stateless; this
adds an indexed tier. Closes the portfolio's no-production-SQL gap with a real relational +
vector schema, and is the first place the M5 `source_version` provenance field gets richly
populated (ChEMBL release, UniProt snapshot — formalized in a corpus manifest).

**Scope:** one target family — the **human kinome** (~500 kinases) and ChEMBL small molecules
with measured bioactivity against them (~tens of thousands of compounds, ~50–100k activity
records). Fit-for-purpose by design; whole store stays under ~1 GB. Full ChEMBL, multiple
families, and per-residue embeddings are explicitly **out of scope for v1**.

### Locked architectural decisions

1. **No `torch` on the request path — all heavy ML is offline.** ESM-2 protein embeddings are
   computed once in the ingestion tier. The serving tier never loads a model. Packaging enforces
   this: ML deps (`torch`, `fair-esm`/`transformers`) live in an **optional `ingest` dependency
   group**, never in the core/serving install.
   - **CPU-only, zero-cost, keyless.** The one-time ~500-protein embedding batch runs **on CPU** —
     a plain `uv run` ingest step, no GPU, no notebook, no cloud. We deliberately do *not* bake in
     a GPU/Colab workflow: GPU access is environment-specific, and a CPU batch over ~500 proteins
     is a reasonable one-time wait that runs anywhere `uv` runs. ESM-2 weights download from a
     public URL (HuggingFace `transformers` / `fair-esm`) with **no API key or account**. The
     whole v0.6.0 path — local Docker Postgres, public model + ChEMBL/UniProt downloads — needs no
     money and no secrets. *(Model size sets the CPU runtime and the schema vector dim — decided at
     M1; the schema currently reserves `vector(1280)` for ESM-2-650M but the corpus is empty, so a
     smaller/faster variant can still be chosen.)*
2. **Chemical similarity = Morgan/Tanimoto (primary), ChemBERTa = stretch + measured baseline.**
   Reuse the existing `tanimoto_similarity` (Morgan radius-2) and InChIKey standardization from
   v0.5.0 — CPU, in-process, no GPU, no request-path model load. This removes the larger GPU
   batch *and* keeps the compound path torch-free. ChemBERTa dense embeddings become a
   **stretch goal whose deliverable is a retrieval-quality comparison vs Tanimoto** (the eval
   artifact), not a request-path dependency. (Consider the RDKit Postgres cartridge `bfp` type
   for in-DB Tanimoto so the hybrid query is one SQL round-trip.)
3. **Protein search over *stored* vectors only.** Serving accepts inputs that resolve to a
   stored embedding (gene symbol / UniProt accession). Embedding a novel raw sequence at
   request time is out of scope for v1 (it would pull ESM-2 onto the request path). This tool
   is framed as a *technique demonstration* — it overlaps existing `get_structural_homologs`
   (Foldseek) on utility.
4. **Exact search at this corpus size; HNSW documented, not reflexive.** ~500 protein vectors
   and tens of thousands of compound vectors are small enough that exact/flat cosine is
   sub-millisecond. Use **SQL-filter-then-exact-kNN** to sidestep the filtered-HNSW recall
   trap, and document the threshold at which HNSW becomes necessary. The vector layer is a
   deliberate *demonstration* of the pattern at portfolio scale, not a performance necessity —
   stated honestly.
5. **ChEMBL: route by task.** Bulk ingestion pulls from the **EBI ChEMBL Postgres/SQLite dump**
   (fast, SQL-sliceable); the existing REST `ChEMBLClient` stays for interactive single-target
   lookups. Right tool for each job.
6. **The corpus DB is optional.** With no `GENESIS_CORPUS_DSN` configured, the server still
   starts and all 38 live-API tools work unchanged; `corpus_*` tools return a typed
   `UpstreamUnavailable` error explaining the corpus isn't configured. Preserves the zero-config
   experience and keeps the existing test suite green.
7. **Same agent contract as everything else.** `corpus_*` tools emit the M5 provenance envelope
   + M6 typed errors; the three confidence signals (vector/Tanimoto similarity, ChEMBL assay
   confidence, `pchembl_value`) are **structured fields in the JSON `data`**, never collapsed
   into one fake score. `readOnlyHint=True`, `openWorldHint=False` (bounded corpus — a
   deliberate contrast with the live-API tools' `openWorldHint=True`).

### Tiers
- **Ingestion (offline)** — `genesis_bio_mcp.ingest` CLI (or `scripts/`), not on the request
  path. Extract (kinase targets+sequences via UniProt/ChEMBL; bioactivities via the ChEMBL
  dump) → Transform (RDKit canonical SMILES + InChIKey + Morgan FP; ESM-2 mean-pooled protein
  vectors, computed once on CPU) → Load (rows + vectors → Postgres; write the `corpus_manifest`).
  Idempotent and re-runnable.
- **Serving (FastMCP)** — read-only `corpus_*` tools query the store at request time via an
  async driver (`asyncpg`), DSN from env. Same Pydantic V2 / `uv` / `ruff` / markdown-out /
  provenance+typed-error conventions.

### Data model
`targets` (target_chembl_id PK, uniprot_accession, gene_symbol, pref_name, organism, sequence,
kinase_group, `sequence_embedding vector(1280)` ESM-2 mean-pooled, source/version/retrieved_at);
`compounds` (molecule_chembl_id PK, canonical_smiles, inchi, inchikey, mol_weight,
**`morgan_fp` ECFP4 — primary similarity**, `chem_embedding vector(768)` ChemBERTa *optional/
stretch*, source/version/retrieved_at); `activities` (activity_id PK, FKs, standard_type/value/
units, pchembl_value, assay_chembl_id, assay_confidence_score, doc_chembl_id, source/version/
retrieved_at); `corpus_manifest` (build timestamp, ChEMBL release, UniProt snapshot, embedding
model+version per modality, row counts, family scope).

### Tools (all read-only, closed-world, paginated, provenance + typed errors)
- **`corpus_describe`** — return the manifest (coverage + provenance + build date). The FAIR/
  transparency entry point.
- **`corpus_search_targets_by_sequence`** — gene/UniProt → stored ESM-2 vector → exact cosine
  kNN over targets; returns top-k targets + per-target bioactivity-coverage count.
- **`corpus_search_compounds`** *(headline hybrid tool)* — SQL filters (target gene/UniProt,
  standard_type, min_pchembl, max_mol_weight) + optional Morgan/Tanimoto similarity from a query
  SMILES; returns compounds with activity values, provenance IDs, and the three confidence
  fields. Paginated.
- **`corpus_find_similar_compounds`** — query SMILES → Morgan/Tanimoto kNN; returns analogs +
  each analog's known target activities in the corpus.

### Testing / CI
Ephemeral Postgres in CI via `testcontainers` / `pytest-postgresql` (honest integration tests,
no live network). Existing respx-mocked suite is untouched because the corpus layer is optional.

### Milestones (one atomic PR each)
- **M0 — Schema + infra.** ✅ *shipped (PR #48):* dockerized pgvector (`docker compose`), schema,
  manifest table, optional async DB pool wired into the lifespan (degrades gracefully when
  unset), `asyncpg`/`pgvector` deps, CI Postgres+pgvector service, `corpus_describe`.
- **M1 — Ingest targets.** ✅ *shipped (PRs #51, #52):* chose **ESM-2-150M (640-dim)** for CPU
  runtime; `corpus_search_targets_by_sequence` (stored-vector cosine kNN) + the CPU kinome
  ingester (UniProt KW-0418 → ESM-2 mean-pool → load). `targets` keyed on UniProt accession.
- **M2 — Ingest compounds + activities.** ✅ *shipped (PRs #53, #54):* `corpus_find_similar_compounds`
  (Morgan/Tanimoto via pgvector bit-Jaccard) + the ChEMBL compound/activity ingester (REST,
  bounded). Activities keyed on UniProt accession + molecule id.
- **M3 — Headline hybrid tool.** ✅ *shipped (PR #55):* `corpus_search_compounds` — SQL filters +
  optional Tanimoto in one query, paginated, three separate confidence signals.
- **M4 — Eval + docs.** ✅ *shipped:* `docs/corpus-eval.md` (10 verifiable QA pairs) + a runnable
  `scripts/corpus_demo.py`; README/docs sections.
- **Stretch.** ChemBERTa embeddings + a Tanimoto-vs-ChemBERTa retrieval-quality comparison.
  (A hosted free-tier demo, e.g. Neon, is *optional* — local Docker Postgres is the default and
  needs no account; only pursue a hosted demo if an always-on live demo is wanted.)

---

## v0.8.0 — Small-molecule depth *(committed; eval-driven)*

The four milestones below are **the four capability gaps the eval harness flagged**
([docs/eval.md](eval.md)) — ADMET/PK, kinome selectivity, patents/FTO, and retrosynthesis.
The eval picked this theme; each milestone exists to flip its probe from GAP → pass.

**Principles (the project's signature):**
- **Keyless / local-compute first** — ship a pure-RDKit tier with no API key, no paid service,
  CPU-only; gate heavier ML behind an *optional* dependency group (never on the serving path).
- **Honest confidence** — every prediction carries an explicit applicability/confidence caveat.
  ADMET, selectivity, and patent calls are *screening signals*, not ground truth — the same
  discipline as the evidence profile and the eval itself.
- **Reuse assets** — RDKit, the ChEMBL client, the kinome corpus (ESM-2 neighbors + activities),
  Tanimoto, and the `_fmt` provenance/typed-error envelope.
- **Each milestone re-runs the eval** to watch its probe flip — the harness is the acceptance test.

### M1 — ADMET / developability panel + synthesizability *(closes `probe-imatinib-admet` + the synthesizability probe)*
- **Tier A (keyless, pure RDKit):** `predict_admet` — ESOL aqueous solubility (Delaney), GI-absorption
  / passive-permeability and BBB heuristics, a **hERG** liability flag (basic amine + high logP), a
  **CYP3A4** metabolic-liability flag, and Brenk/PAINS structural alerts — each with an honest "rule-based
  alert, not a validated prediction" caveat. **SAscore** (ships with RDKit Contrib — confirmed) for
  synthetic accessibility.
- **Tier B (optional `admet` dep group, later):** quantitative ML endpoints (hERG, CYP1A2/2C9/2D6/3A4,
  clearance, t½) via ADMET-AI (Chemprop) or sklearn-on-TDC — returned with applicability-domain notes.
- New `tools/admet.py` (pure compute) + `ADMETProfile` model + `predict_admet` tool.

### M2 — Kinome selectivity / off-target *(closes `probe-panraf-selectivity`)*
The **corpus payoff**. `assess_selectivity` returns **measured** off-targets (ChEMBL + corpus
`activities`) and **predicted** off-target risk (ESM-2 sequence-neighbor kinases for a target;
Tanimoto-to-known-binders, SEA-style, for a compound) — honestly separated, with a coverage note.
Measured layer is keyless; the predicted layer wants a built corpus (gated `requires: ["corpus"]`).

### M3 — Patents / freedom-to-operate *(closes `probe-krasg12c-patents`)*
`search_patents` over **SureChEMBL** (compound→patent, chemistry-native) + optional PatentsView/USPTO
landscape, carrying a **permanent "informational, not FTO legal advice" disclaimer** (the OpenFDA
disclaimer pattern). Scoped honestly: chemistry-in-patents, US-centric. Highest external-data friction →
sequenced last of the keyless-ish set.

### M4 (stretch) — Retrosynthetic route planning *(closes the route-planning probe)*
AiZynthFinder (AstraZeneca, MIT) in an optional dep group with downloaded templates. Heavy, lowest ROI
→ explicit stretch after M1–M3.

**Sequence:** M1 → M2 → M3, then the M4 stretch. After each, re-run `python -m genesis_bio_mcp.evals`.

## v0.7.0 and beyond *(directional, not committed)*

- **Structure-based methods** — pocket detection / druggability (fpocket-style), lightweight
  docking or binding-affinity scoring.
- **Generative / design** — scaffold hopping, R-group enumeration.
- **Proprietary-data seam (`Architecture4Insight`)** — `docs/private-data-sources.md` + an
  auth scaffold (API-key/OAuth per client) and config-driven enable/disable, so an org can
  plug an internal compound registry / assay DB / ELN-LIMS behind the *same* MCP interface and
  tool contract. (The v0.6.0 corpus layer is the natural substrate for this.)
- **Indication→tissue mapping (EFO→UBERON)**, Reactome parent-pathway labeling, and the
  remaining deferred v0.3.x items.

---

_Roadmap is a living document — milestones may be re-scoped as upstream APIs and priorities evolve._
