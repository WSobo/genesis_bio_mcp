"""Task-based evaluation harness for genesis-bio-mcp.

A small, curated set of drug-discovery questions with *verifiable* expectations,
scored in two layers:

- **Deterministic (keyless):** call each item's declared tools directly and grade
  the output — proves whether the answer is *retrievable* with the current tools.
- **Agentic (optional, needs ANTHROPIC_API_KEY):** run the workflow agent end-to-end
  and score tool-selection + answer correctness.

The report's **Capability gaps** section — questions no current tool can answer —
is the evidence-driven roadmap for what to build next.

Run it:  ``uv run python -m genesis_bio_mcp.evals --out report.md [--agentic]``
"""

from genesis_bio_mcp.evals.dataset import EvalItem, load_dataset
from genesis_bio_mcp.evals.graders import GradeResult, grade
from genesis_bio_mcp.evals.report import EvalReport, ItemResult
from genesis_bio_mcp.evals.runner import run_eval

__all__ = [
    "EvalItem",
    "load_dataset",
    "grade",
    "GradeResult",
    "EvalReport",
    "ItemResult",
    "run_eval",
]
