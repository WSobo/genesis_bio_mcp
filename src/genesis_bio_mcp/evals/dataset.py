"""Eval dataset: typed items + a validating loader.

The dataset (``dataset.json``) is the machine source of truth for the harness.
Each item declares the deterministic ``tool_calls`` that should surface the answer,
the ``expected_tools`` an agent ought to choose, and a ``grader`` (see graders.py).
"""

from __future__ import annotations

import json
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field

_DATASET_PATH = Path(__file__).parent / "dataset.json"


class ToolCall(BaseModel):
    """One deterministic tool invocation: a registry tool name + its kwargs."""

    model_config = ConfigDict(extra="forbid")
    name: str
    args: dict = Field(default_factory=dict)


class EvalItem(BaseModel):
    """A single evaluation item with a verifiable expectation."""

    model_config = ConfigDict(extra="forbid")
    id: str
    question: str = Field(description="Natural-language question (drives the agentic layer)")
    category: str = Field(description="Tool category for coverage stats")
    tool_calls: list[ToolCall] = Field(
        default_factory=list,
        description="Deterministic path: tools to call (in order) to surface the answer",
    )
    expected_tools: list[str] = Field(
        default_factory=list,
        description="Tools an agent should select (tool-selection ground truth)",
    )
    grader: dict = Field(description="Grader spec applied to tool/agent output (see graders.py)")
    requires: list[str] = Field(
        default_factory=list, description="Gating tags, e.g. 'corpus', 'anthropic_api_key'"
    )
    probe: bool = Field(
        default=False,
        description="True = no current tool can answer this; a deterministic fail is a CAPABILITY GAP, not a regression",
    )

    # --- Optional authoring metadata (see the mcp-eval-authoring skill) ------------------------
    # All default-empty, so existing items validate unchanged. None of these affect scoring; they
    # record how an item was authored/checked and make the suite legible to reviewers and to
    # future re-baselining. `EvalItem` is extra="forbid", so these must be declared to be usable.
    notes: str = Field(
        default="",
        description="Author notes: the correct path, the naive path a hurried agent would take, and the tolerance rationale",
    )
    source_version: dict = Field(
        default_factory=dict,
        description="Upstream versions observed at authoring, e.g. {'chembl': '35', 'observed_at': '2026-07-01'} — lets a later failure be diagnosed as drift vs regression",
    )
    volatile: bool = Field(
        default=False,
        description="True if the grader targets a volatile upstream scalar; expect to re-baseline on upstream releases",
    )
    baseline_expectation: str | None = Field(
        default=None,
        description="Expected tools-off outcome ('pass' | 'fail' | 'fail (stochastic)'); a 'pass' flags a shortcut to rewrite (measured via `--baseline`)",
    )
    canary: str = Field(
        default="",
        description="Canary string so the item is detectable if it leaks into a future training corpus",
    )
    replicates: int = Field(
        default=1, ge=1, description="K — recommended replicate count for this item (Gate 4)"
    )
    failure_modes: list[str] = Field(
        default_factory=list,
        description="Observed failure-mode tags, e.g. 'tool_not_found', 'output_misread' (Gate 5)",
    )
    failure_mode_taxonomy_version: str = Field(
        default="", description="Version of the failure-mode vocabulary these tags use"
    )


def load_dataset(path: Path | None = None) -> list[EvalItem]:
    """Load + validate the dataset, enforcing unique ids."""
    raw = json.loads((path or _DATASET_PATH).read_text())
    items = [EvalItem.model_validate(x) for x in raw]
    ids = [i.id for i in items]
    dupes = sorted({x for x in ids if ids.count(x) > 1})
    if dupes:
        raise ValueError(f"duplicate eval ids: {dupes}")
    return items
