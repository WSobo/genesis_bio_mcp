"""Tests for the eval harness — all keyless / offline, so they run in CI.

The live eval run (``python -m genesis_bio_mcp.evals``) hits real APIs and is
exercised on demand; here we test the machinery: graders, dataset validity,
runner internals, report rendering, and the agent trace capture (mocked client).
"""

from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import AsyncMock, MagicMock

from genesis_bio_mcp.evals.dataset import EvalItem, ToolCall, load_dataset
from genesis_bio_mcp.evals.graders import grade
from genesis_bio_mcp.evals.report import EvalReport, ItemResult
from genesis_bio_mcp.evals.runner import _availability, _run_tool_calls
from genesis_bio_mcp.workflow_agent import ToolSpec, build_tool_registry, run_agent_loop_traced
from tests.test_workflow_agent import _make_response, _text_block, _tool_use_block


def _spec(name: str, fn, props: dict) -> ToolSpec:
    return ToolSpec(
        name=name,
        description="d",
        input_schema={"type": "object", "properties": props},
        tool_category="x",
        use_when="w",
        fn=fn,
    )


# ---------------------------------------------------------------------------
# Graders
# ---------------------------------------------------------------------------


def test_grader_contains_case_insensitive():
    assert grade("UniProt P04626 here", {"type": "contains", "value": "p04626"}).passed
    assert not grade("nothing", {"type": "contains", "value": "P04626"}).passed


def test_grader_all_of():
    assert grade("BRAF and KRAS", {"type": "all_of", "value": ["BRAF", "KRAS"]}).passed
    assert not grade("BRAF only", {"type": "all_of", "value": ["BRAF", "KRAS"]}).passed


def test_grader_regex():
    assert grade("MAPK cascade", {"type": "regex", "value": "MAPK|RAF"}).passed
    assert not grade("nothing here", {"type": "regex", "value": "MAPK|RAF"}).passed


def test_grader_numeric_value_tolerance():
    assert grade("MW 180.16 g/mol", {"type": "numeric", "value": 180.0, "tolerance": 0.5}).passed
    assert not grade("MW 200", {"type": "numeric", "value": 180.0, "tolerance": 0.5}).passed


def test_grader_numeric_min_with_near_label():
    assert grade(
        "best pChEMBL 9.5 reported", {"type": "numeric", "near": "pChEMBL", "min": 8.0}
    ).passed
    assert not grade("best pChEMBL 6.0", {"type": "numeric", "near": "pChEMBL", "min": 8.0}).passed
    # label absent → fail even though a qualifying number exists elsewhere
    assert not grade("9.9 somewhere", {"type": "numeric", "near": "pChEMBL", "min": 8.0}).passed


def test_grader_set_recall_at_k():
    canon = ["ABL1", "LYN", "FYN"]
    # default recall_at_k = 1.0 → all must be present (generalizes all_of)
    assert grade("kinases ABL1, Lyn, and FYN", {"type": "set", "canonical": canon}).passed
    assert not grade("ABL1 and LYN only", {"type": "set", "canonical": canon}).passed
    # threshold 0.66 → 2 of 3 passes, 1 of 3 fails
    assert grade("ABL1 and LYN", {"type": "set", "canonical": canon, "recall_at_k": 0.66}).passed
    assert not grade("just ABL1", {"type": "set", "canonical": canon, "recall_at_k": 0.66}).passed
    # case-insensitive; empty canonical fails rather than crashing
    assert grade("abl1 lyn fyn", {"type": "set", "canonical": canon}).passed
    assert not grade("anything", {"type": "set", "canonical": []}).passed
    # detail reports recall fraction + what's missing
    r = grade("ABL1 only", {"type": "set", "canonical": canon, "recall_at_k": 0.66})
    assert "1/3" in r.detail and "LYN" in r.detail


def test_grader_unknown_type_and_none_output():
    assert not grade("x", {"type": "bogus"}).passed
    assert not grade(None, {"type": "contains", "value": "x"}).passed


# ---------------------------------------------------------------------------
# Dataset validity (the real shipped dataset)
# ---------------------------------------------------------------------------


def test_dataset_loads_with_unique_ids():
    items = load_dataset()
    assert len(items) >= 20
    ids = [i.id for i in items]
    assert len(ids) == len(set(ids))


def test_dataset_tool_calls_valid_against_registry():
    """Every tool_call name + every arg key must exist in the live registry schema."""
    reg = build_tool_registry(SimpleNamespace())
    for item in load_dataset():
        for tc in item.tool_calls:
            assert tc.name in reg, f"{item.id}: unknown tool {tc.name}"
            props = set((reg[tc.name].input_schema or {}).get("properties", {}))
            bad = [k for k in tc.args if k not in props]
            assert not bad, f"{item.id}: {tc.name} bad args {bad}"
        for et in item.expected_tools:
            assert et in reg, f"{item.id}: unknown expected_tool {et}"


def test_probe_items_have_no_tool_calls():
    for item in load_dataset():
        if item.probe:
            assert item.tool_calls == [], f"{item.id}: a probe must have no tool_calls"


def test_eval_item_authoring_metadata_defaults_and_roundtrip():
    """The optional mcp-eval-authoring metadata fields default empty (so pre-existing items
    validate unchanged) and accept values when provided."""
    bare = EvalItem(id="a", question="q", category="c", grader={"type": "contains", "value": "x"})
    assert bare.notes == "" and bare.canary == ""
    assert bare.source_version == {} and bare.failure_modes == []
    assert bare.volatile is False and bare.baseline_expectation is None
    assert bare.replicates == 1 and bare.failure_mode_taxonomy_version == ""

    rich = EvalItem(
        id="b",
        question="q",
        category="druggability",
        grader={"type": "numeric", "near": "QED", "value": 0.46, "tolerance": 0.03},
        notes="correct: RDKit QED; naive: guess ~0.5",
        source_version={"rdkit": "2024.03", "observed_at": "2026-07-23"},
        baseline_expectation="fail (stochastic)",
        canary="canary GUID 00000000-0000-0000-0000-000000000000",
        replicates=3,
        failure_modes=["output_misread"],
        failure_mode_taxonomy_version="v1",
    )
    assert rich.replicates == 3
    assert rich.baseline_expectation == "fail (stochastic)"
    assert rich.model_dump()["failure_modes"] == ["output_misread"]


def test_eval_item_still_rejects_unknown_field():
    """extra='forbid' must still reject a genuinely unknown field (typo guard)."""
    import pytest
    from pydantic import ValidationError

    with pytest.raises(ValidationError):
        EvalItem(
            id="c",
            question="q",
            category="c",
            grader={"type": "contains", "value": "x"},
            notez="typo",  # not a real field
        )


# ---------------------------------------------------------------------------
# Runner internals (no network)
# ---------------------------------------------------------------------------


async def test_run_tool_calls_concatenates_outputs():
    async def fake(gene_symbol):
        return f"result for {gene_symbol}"

    reg = {"t": _spec("t", fake, {"gene_symbol": {"type": "string"}})}
    item = EvalItem(
        id="x",
        question="q",
        category="c",
        tool_calls=[ToolCall(name="t", args={"gene_symbol": "BRAF"})],
        grader={"type": "contains", "value": "BRAF"},
    )
    out, err = await _run_tool_calls(item, reg)
    assert "BRAF" in out
    assert err == ""


async def test_run_tool_calls_unknown_tool_reports_error():
    item = EvalItem(
        id="y",
        question="q",
        category="c",
        tool_calls=[ToolCall(name="missing", args={})],
        grader={"type": "contains", "value": "z"},
    )
    out, err = await _run_tool_calls(item, {})
    assert "unknown tool" in err


def test_availability_gates_on_requires():
    item = EvalItem(
        id="x",
        question="q",
        category="c",
        grader={"type": "contains", "value": "z"},
        requires=["corpus"],
    )
    ok, reason = _availability(item, has_key=True, has_corpus=False)
    assert not ok and "corpus" in reason
    ok2, _ = _availability(item, has_key=True, has_corpus=True)
    assert ok2


# ---------------------------------------------------------------------------
# Report rendering
# ---------------------------------------------------------------------------


def test_report_markdown_surfaces_gaps_and_summary():
    results = [
        ItemResult(id="ok", category="gene_annotation", probe=False, question="q", det_passed=True),
        ItemResult(
            id="gap", category="druggability", probe=True, question="ADMET?", det_passed=False
        ),
    ]
    rep = EvalReport(generated_at="2026-01-01T00:00:00Z", results=results, agentic=False)
    md = rep.to_markdown()
    assert "Capability gaps" in md
    assert "ADMET?" in md
    s = rep.to_json_dict()["summary"]
    assert s["deterministic_passed"] == 1
    assert s["deterministic_total"] == 1
    assert s["capability_gaps"] == 1


# ---------------------------------------------------------------------------
# Agent trace capture (mocked Anthropic client)
# ---------------------------------------------------------------------------


def _mock_client(responses: list) -> MagicMock:
    client = MagicMock()
    client.messages = MagicMock()
    client.messages.create = AsyncMock(side_effect=responses)
    return client


async def test_traced_loop_records_tool_calls_and_terminates():
    async def fake(gene_name):
        return "ERBB2 P04626"

    reg = {"resolve_gene": _spec("resolve_gene", fake, {"gene_name": {"type": "string"}})}
    step1 = _make_response(
        "tool_use", [_tool_use_block("resolve_gene", {"gene_name": "HER2"}, "t1")]
    )
    final = _make_response("end_turn", [_text_block("The accession is P04626.")])
    run = await run_agent_loop_traced("q", reg, client=_mock_client([step1, final]))

    assert run.tool_calls == ["resolve_gene"]
    assert run.stop_reason == "end_turn"
    assert run.iterations == 2
    assert "P04626" in run.final_text


async def test_traced_loop_caps_at_max_iterations():
    async def fake(gene_name):
        return "x"

    reg = {"resolve_gene": _spec("resolve_gene", fake, {"gene_name": {"type": "string"}})}
    step = _make_response(
        "tool_use", [_tool_use_block("resolve_gene", {"gene_name": "HER2"}, "t1")]
    )
    client = MagicMock()
    client.messages = MagicMock()
    client.messages.create = AsyncMock(return_value=step)  # never ends
    run = await run_agent_loop_traced("q", reg, client=client, max_iterations=3)

    assert run.stop_reason == "max_iterations"
    assert run.iterations == 3
    assert run.tool_calls == ["resolve_gene", "resolve_gene", "resolve_gene"]


# ---------------------------------------------------------------------------
# Tools-off (prior-only) baseline + uplift
# ---------------------------------------------------------------------------


async def test_traced_loop_empty_registry_is_prior_only_no_tools():
    """The tools-off baseline: an empty registry sends NO `tools` param (the API rejects
    `tools=[]`) and the model answers from memory in a single turn."""
    final = _make_response("end_turn", [_text_block("The symbol is PTGS2.")])
    client = _mock_client([final])
    run = await run_agent_loop_traced("q", {}, client=client)

    assert run.tool_calls == []
    assert run.stop_reason == "end_turn"
    assert "PTGS2" in run.final_text
    assert "tools" not in client.messages.create.call_args.kwargs


def test_report_uplift_surfaces_tool_value():
    results = [
        ItemResult(
            id="tool-dep",
            category="druggability",
            probe=False,
            question="q",
            det_passed=True,
            agentic_ran=True,
            agentic_passed=True,
            baseline_ran=True,
            baseline_passed=False,  # only answerable WITH tools
        ),
        ItemResult(
            id="shortcut",
            category="gene_annotation",
            probe=False,
            question="q",
            det_passed=True,
            agentic_ran=True,
            agentic_passed=True,
            baseline_ran=True,
            baseline_passed=True,  # answerable WITHOUT tools
        ),
    ]
    rep = EvalReport(
        generated_at="2026-01-01T00:00:00Z", results=results, agentic=True, baseline=True
    )

    s = rep.to_json_dict()["summary"]
    assert s["tools_on_passed"] == 2
    assert s["tools_off_passed"] == 1
    assert s["uplift_net"] == 1
    assert s["tool_dependent"] == 1
    assert s["shortcuts"] == 1

    md = rep.to_markdown()
    assert "Uplift" in md
    assert "Tool-value analysis" in md
    assert "tool-dep" in md  # only-answerable-with-tools
    assert "shortcut" in md  # a shortcut to rewrite/cut
    assert "Tools-off" in md  # per-item results column

    row = next(r for r in rep.to_json_dict()["results"] if r["id"] == "shortcut")
    assert row["baseline_ran"] and row["baseline_passed"]


def test_report_results_table_columns_align():
    """Regression: each results row must have exactly as many cells as the header. A prior
    ``row[:-1]`` slice dropped a column separator and merged the latency/agentic cells."""
    r = ItemResult(
        id="x",
        category="c",
        probe=False,
        question="q",
        det_passed=True,
        det_latency_s=1.0,
        agentic_ran=True,
        agentic_passed=True,
        tool_selection_ok=True,
        baseline_ran=True,
        baseline_passed=False,
    )
    for baseline in (False, True):
        md = EvalReport(
            generated_at="t", results=[r], agentic=True, baseline=baseline
        ).to_markdown()
        table = md.split("## Results")[1].strip().splitlines()
        header, _sep, body = table[0], table[1], table[2]
        assert body.count("|") == header.count("|") == _sep.count("|")
