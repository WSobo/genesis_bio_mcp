# Evaluation harness

A task-based eval that measures whether the MCP's tools (and the workflow agent that
selects them) actually answer realistic drug-discovery questions — and, just as
importantly, **surfaces what they can't answer yet**. It replaces ad-hoc smoke-testing
with an evidence-driven loop: the report's **Capability gaps** section is the roadmap
for what to build next.

The dataset lives at [`src/genesis_bio_mcp/evals/dataset.json`](../src/genesis_bio_mcp/evals/dataset.json)
— ~20 questions with *verifiable* expectations plus a handful of deliberate
capability **probes** (questions no current tool can answer).

Every non-probe item is authored to be **tool-necessary**: the answer hinges on a precise
value or identifier the tool returns and a strong model cannot recall unaided (an InChIKey,
an InterPro/MONDO/Reactome accession, a gnomAD constraint value, an AlphaFold pLDDT, a
computed QED/TPSA/SAscore), and the target value never appears in the question. This is
deliberate: a question a model answers cold (e.g. "HER2's UniProt accession → P04626")
scores a pass forever and measures training data, not the server. The harness's tools-off
`--baseline` mode exists to *detect* such shortcut items automatically — an item a no-tools
run still passes is a broken item, not a hard one.

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

**3. Uplift baseline (`--baseline`, implies `--agentic`).** Runs each item a *second* time
through the same agent with an **empty tool registry** — a *prior-only* run that must answer
from the model's own knowledge. Grading both runs with the same grader gives the headline
**uplift = tools-on − tools-off**, and the report's *Tool-value analysis* section splits items
into **only-answerable-with-tools** vs. **shortcuts** (the model already knew the answer, so
they measure prior knowledge, not tool value — rewrite or cut them). This is the check a
deterministic `N/N` score cannot give you: a passing item only proves the server *helps* if a
tools-off agent *fails* it. (The deterministic layer runs the declared tools then greps their
own output, so it always passes — it verifies the data path, not tool value.)

## Running it

```bash
# Deterministic only (no key needed). Prints a Markdown report; --out also writes JSON.
uv run python -m genesis_bio_mcp.evals --out eval_report.md

# Add the agentic layer (requires ANTHROPIC_API_KEY in the environment)
uv run python -m genesis_bio_mcp.evals --agentic --out eval_report.md

# Add the tools-off uplift baseline — the headline "do the tools help?" number (implies --agentic)
uv run python -m genesis_bio_mcp.evals --baseline --out eval_report.md

# Replicate the LLM layers K times for confidence intervals (one run is an anecdote)
uv run python -m genesis_bio_mcp.evals --baseline --replicates 3 --out eval_report.md

# Run a single category
uv run python -m genesis_bio_mcp.evals --filter druggability
```

`--replicates K` re-runs the **agentic and baseline layers K times** per item (the deterministic
layer is deterministic, so it runs once). Each item's pass **rate** (`passes/K`) is aggregated
across items into a **t-based 95% CI** — the *item* is the sampling unit, so the interval reflects
both LLM sampling variance and item-to-item spread. The report then shows `agentic rate`,
`tool-selection`, `tools-off baseline`, and `uplift (Δ)` as `mean [lo, hi]` instead of a single
point estimate. This is what makes an uplift number defensible rather than anecdotal — the
tools-off residual is genuinely stochastic run to run, so a single pass is noise.

## Reading the report

- **Summary** — deterministic pass rate, how many probes the current tools could answer
  (the rest are gaps), and (agentic) tool-selection accuracy.
- **Coverage by category** — which tool categories are exercised and how they score.
- **Capability gaps — what to build next** — the headline. Failed probes are *expected*
  gaps; any non-probe failure is flagged separately to investigate (a miscalibrated grader
  or a real data-layer regression).
- **Tool-value analysis (uplift)** — only with `--baseline`. Lists items only answerable with
  tools (genuine value) and shortcuts a tools-off agent already passes (rewrite or cut). If this
  section is mostly shortcuts, the dataset is measuring the model, not the server.

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
optionally `near` a label), and `set` — recall over a `canonical` list that must meet
`recall_at_k` (default 1.0), the right grader for list answers (off-target kinases, pathway/ID
sets) and a soft generalization of `all_of`. A `"probe": true` item has no `tool_calls` — a
deterministic failure on it is a capability gap, not a regression. Add `"requires": ["corpus"]`
to gate an item on a built corpus.

The harness machinery (graders, dataset validity, report rendering, agent-trace capture)
is unit-tested in [`tests/test_evals.py`](../tests/test_evals.py) — keyless and CI-safe.

## Tool-necessity, measured

The dataset was rebuilt (2026-07) so every non-probe item is **tool-necessary** — its answer
is a precise value/identifier the tool returns and a strong model cannot recall unaided (see
the design note above). We don't just assert this; we measure it with the tools-off baseline:

| | Deterministic (tools) | Tools-off (prior-only) |
|---|---:|---:|
| Rebuilt dataset (18 items) | 18/18 | ~3/18 |
| *Previous dataset (for contrast)* | *27/27* | *25/27* |

The previous suite was ~93% answerable with the tools unplugged — it graded facts a strong
model already knows (`P04626`, `PTGS2`) or matched tool-format headers, so its high pass rate
measured **training data, not the server**. The rebuilt suite inverts that: pull the tools and
the score collapses to ~17%, so a pass now largely reflects the tool *adding* something.
Authoring items to hit this bar — no-tool baseline, discriminating graders, drift-stable
answers — follows the `mcp-eval-authoring` methodology.

The residual ~3/18 is **stochastic** (which items leak varies run to run): against a saturated
model a few answers are still estimable within tolerance (a clean molecule's QED to ±0.03) or
mentionable (a related accession). Closing that last gap needs replicated baselines and
structured-answer graders (precision/recall on a returned set, an exact typed field) rather
than the current substring/numeric graders — the natural next increment to the harness.

Two deliberate **capability probes** remain unmet — patents/FTO and retrosynthesis — the same
gaps the roadmap tracks as v0.8.0 M3/M4.

For the corpus-specific eval set, see [corpus-eval.md](corpus-eval.md).
