"""Eval runner: builds a live client graph and scores items in two layers.

Layer 1 (deterministic, keyless): call each item's declared ``tool_calls`` directly
through the registry and grade the concatenated output. Layer 2 (agentic, optional,
needs ANTHROPIC_API_KEY): run the workflow agent and score tool-selection + answer.

Makes live HTTP calls — an on-demand harness, not a CI job.
"""

from __future__ import annotations

import os
import time

import httpx

from genesis_bio_mcp.config.settings import settings
from genesis_bio_mcp.evals.dataset import EvalItem
from genesis_bio_mcp.evals.graders import grade
from genesis_bio_mcp.evals.report import EvalReport, ItemResult
from genesis_bio_mcp.server import build_app_state
from genesis_bio_mcp.workflow_agent import build_tool_registry, run_agent_loop_traced

_EVAL_HEADERS = {
    "User-Agent": "genesis-bio-mcp-eval/0.1 (+https://github.com/WSobo/genesis-bio-mcp)"
}


def _availability(item: EvalItem, *, has_key: bool, has_corpus: bool) -> tuple[bool, str]:
    """Decide whether an item can run given the environment; returns (ok, skip_reason)."""
    for req in item.requires:
        if req == "anthropic_api_key" and not has_key:
            return False, "no ANTHROPIC_API_KEY"
        if req == "corpus" and not has_corpus:
            return False, "no corpus DSN"
    return True, ""


async def _run_tool_calls(item: EvalItem, registry: dict) -> tuple[str, str]:
    """Execute an item's tool_calls in order; return (concatenated_output, error)."""
    outputs: list[str] = []
    for tc in item.tool_calls:
        spec = registry.get(tc.name)
        if spec is None:
            return "\n".join(outputs), f"unknown tool '{tc.name}'"
        try:
            out = await spec.fn(**tc.args)
        except Exception as exc:
            return "\n".join(outputs), f"{tc.name}: {exc}"
        outputs.append(out or "")
    return "\n".join(outputs), ""


async def run_eval(
    items: list[EvalItem],
    *,
    agentic: bool = False,
    max_iterations: int = 8,
    timestamp: str,
) -> EvalReport:
    """Run the eval over ``items`` and return an :class:`EvalReport`."""
    has_key = bool(os.environ.get("ANTHROPIC_API_KEY"))
    results: list[ItemResult] = []

    async with httpx.AsyncClient(
        headers=_EVAL_HEADERS, timeout=settings.httpx_timeout, follow_redirects=True
    ) as client:
        state = build_app_state(
            client
        )  # no DepMap cache → DepMap uses its OT proxy (fine for eval)
        registry = build_tool_registry(state)
        has_corpus = getattr(state, "corpus_pool", None) is not None

        for item in items:
            r = ItemResult(
                id=item.id, category=item.category, probe=item.probe, question=item.question
            )
            ok, reason = _availability(item, has_key=has_key, has_corpus=has_corpus)
            if not ok:
                r.skipped, r.skip_reason = True, reason
                results.append(r)
                continue

            # Layer 1 — deterministic
            t0 = time.perf_counter()
            output, err = await _run_tool_calls(item, registry)
            r.det_latency_s = time.perf_counter() - t0
            r.det_error = err
            g = grade(output, item.grader)
            r.det_passed = g.passed
            r.det_detail = err or g.detail

            # Layer 2 — agentic (optional)
            if agentic and has_key:
                run = await run_agent_loop_traced(
                    item.question, registry, max_iterations=max_iterations
                )
                r.agentic_ran = True
                r.called_tools = run.tool_calls
                ag = grade(run.final_text, item.grader)
                r.agentic_passed = ag.passed
                r.agentic_detail = ag.detail
                if item.expected_tools:
                    r.tool_selection_ok = all(t in run.tool_calls for t in item.expected_tools)

            results.append(r)

    return EvalReport(generated_at=timestamp, results=results, agentic=agentic)
