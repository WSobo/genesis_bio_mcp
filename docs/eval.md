# Evaluation harness

A task-based eval that measures whether the MCP's tools (and the workflow agent that
selects them) actually answer realistic drug-discovery questions — and, just as
importantly, **surfaces what they can't answer yet**. It replaces ad-hoc smoke-testing
with an evidence-driven loop: the report's **Capability gaps** section is the roadmap
for what to build next.

The dataset lives at [`src/genesis_bio_mcp/evals/dataset.json`](../src/genesis_bio_mcp/evals/dataset.json)
— ~24 questions with *verifiable* expectations plus a handful of deliberate
capability **probes** (questions no current tool can answer).

## Two scoring layers

**1. Deterministic (keyless).** For each item, the runner calls the item's declared
`tool_calls` directly through the tool registry and grades the concatenated output. This
answers *"is the right answer retrievable with the current tools?"* — it needs no API key
and makes live HTTP calls, so it runs on demand (not in CI).

**2. Agentic (optional, needs `ANTHROPIC_API_KEY`).** Runs the full `run_biology_workflow`
agent on the natural-language question and scores two things: **tool selection**
(`expected_tools` ⊆ the tools the agent actually called, via the `AgentRun` trace) and
**answer correctness** (the same grader applied to the agent's final text). Skipped — not
failed — when no key is set.

## Running it

```bash
# Deterministic only (no key needed). Prints a Markdown report; --out also writes JSON.
uv run python -m genesis_bio_mcp.evals --out eval_report.md

# Add the agentic layer (requires ANTHROPIC_API_KEY in the environment)
uv run python -m genesis_bio_mcp.evals --agentic --out eval_report.md

# Run a single category
uv run python -m genesis_bio_mcp.evals --filter druggability
```

## Reading the report

- **Summary** — deterministic pass rate, how many probes the current tools could answer
  (the rest are gaps), and (agentic) tool-selection accuracy.
- **Coverage by category** — which tool categories are exercised and how they score.
- **Capability gaps — what to build next** — the headline. Failed probes are *expected*
  gaps; any non-probe failure is flagged separately to investigate (a miscalibrated grader
  or a real data-layer regression).

## Dataset format

Each item declares the deterministic path, the agent's expected tool selection, and a
grader (see [`evals/graders.py`](../src/genesis_bio_mcp/evals/graders.py)):

```json
{
  "id": "her2-uniprot-accession",
  "question": "HER2 is the clinical name ... give its UniProt accession.",
  "category": "gene_annotation",
  "tool_calls": [
    {"name": "resolve_gene", "args": {"gene_name": "HER2"}},
    {"name": "get_protein_info", "args": {"gene_symbol": "ERBB2"}}
  ],
  "expected_tools": ["resolve_gene", "get_protein_info"],
  "grader": {"type": "contains", "value": "P04626"}
}
```

Grader types: `contains`, `all_of`, `regex`, `numeric` (`value`+`tolerance`, or `min`/`max`,
optionally `near` a label). A `"probe": true` item has no `tool_calls` — a deterministic
failure on it is a capability gap, not a regression. Add `"requires": ["corpus"]` to gate
an item on a built corpus.

The harness machinery (graders, dataset validity, report rendering, agent-trace capture)
is unit-tested in [`tests/test_evals.py`](../tests/test_evals.py) — keyless and CI-safe.

## First run (baseline)

The initial deterministic run scored **23/24** answerable tasks across all eight tool
categories and flagged **4 capability gaps** — none of which any current tool can address:

| Gap | What it needs |
|---|---|
| ADMET / PK prediction (hERG, CYP3A4) | a computed cheminformatics ADMET tool |
| Kinome-wide selectivity profile | off-target / polypharmacology analysis |
| Patents / freedom-to-operate | a competitive-landscape data source |
| Retrosynthesis / synthesizability | a synthesis-planning method |

That table is the evidence behind the next roadmap decision — small-molecule depth leads.

For the corpus-specific eval set, see [corpus-eval.md](corpus-eval.md).
