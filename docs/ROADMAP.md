# Roadmap

Current release: **v0.5.0** — 38 tools across 24 public biomedical data sources (plus local
RDKit cheminformatics and the UMA-Inverse model service), with `prioritize_target` /
`compare_targets` orchestration, the `run_biology_workflow` agent, and the agent-contract
provenance/typed-error envelope on all JSON output.

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

## v0.6.0 and beyond *(directional, not committed)*

- **Computed ADMET (`assess_admet`)** — predicted solubility, permeability, CYP/hERG/tox via an
  open model (e.g. ADMET-AI). Deferred from v0.5.0 because it adds an ML dependency + latency.
- **Structure-based methods** — pocket detection / druggability (fpocket-style), lightweight
  docking or binding-affinity scoring.
- **Generative / design** — scaffold hopping, R-group enumeration.
- **Proprietary-data seam (`Architecture4Insight`)** — `docs/private-data-sources.md` + an
  auth scaffold (API-key/OAuth per client) and config-driven enable/disable, so an org can
  plug an internal compound registry / assay DB / ELN-LIMS behind the *same* MCP interface and
  tool contract. (Enables enterprise adoption without shipping any proprietary data.)
- **Indication→tissue mapping (EFO→UBERON)**, Reactome parent-pathway labeling, and the
  remaining deferred v0.3.x items.

---

_Roadmap is a living document — milestones may be re-scoped as upstream APIs and priorities evolve._
