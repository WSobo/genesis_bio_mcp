"""CLI entrypoint for the eval harness.

    uv run python -m genesis_bio_mcp.evals [--agentic] [--baseline] [--replicates K] [--out report.md] [--filter CATEGORY]

Always runs the deterministic (keyless) layer; ``--agentic`` adds the workflow-agent
layer when ANTHROPIC_API_KEY is set. ``--baseline`` (implies ``--agentic``) additionally
runs each item with NO tools — a prior-only baseline — and reports the **uplift**
(tools-on minus tools-off): the headline "do the tools actually help?" number. ``--replicates K``
re-runs the LLM layers K times per item (deterministic runs once) and reports a t-based 95%
confidence interval over per-item rates — one run is an anecdote. Prints the Markdown report;
``--out`` also writes it plus a sibling ``.json``.
"""

from __future__ import annotations

import argparse
import asyncio
import json
from datetime import UTC, datetime
from pathlib import Path

from genesis_bio_mcp.evals.dataset import load_dataset
from genesis_bio_mcp.evals.runner import run_eval


def main() -> None:
    ap = argparse.ArgumentParser(prog="genesis_bio_mcp.evals", description=__doc__)
    ap.add_argument("--agentic", action="store_true", help="also run the workflow-agent layer")
    ap.add_argument(
        "--baseline",
        action="store_true",
        help="also run a tools-off prior-only baseline and report uplift (implies --agentic)",
    )
    ap.add_argument(
        "--replicates",
        type=int,
        default=1,
        help="re-run the LLM layers K times per item for a confidence interval (agentic only)",
    )
    ap.add_argument("--out", type=str, default=None, help="write report.md (+ .json) to this path")
    ap.add_argument("--filter", type=str, default=None, help="only run items in this category")
    ap.add_argument("--max-iterations", type=int, default=8, help="agent loop cap (agentic only)")
    args = ap.parse_args()

    items = load_dataset()
    if args.filter:
        items = [i for i in items if i.category == args.filter]
        if not items:
            raise SystemExit(f"no eval items in category {args.filter!r}")

    timestamp = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    agentic = args.agentic or args.baseline  # baseline is meaningless without the tools-on run
    report = asyncio.run(
        run_eval(
            items,
            agentic=agentic,
            baseline=args.baseline,
            replicates=args.replicates,
            max_iterations=args.max_iterations,
            timestamp=timestamp,
        )
    )

    md = report.to_markdown()
    print(md)
    if args.out:
        out = Path(args.out)
        out.write_text(md)
        out.with_suffix(".json").write_text(json.dumps(report.to_json_dict(), indent=2))
        print(f"\n[wrote {out} and {out.with_suffix('.json')}]")


if __name__ == "__main__":
    main()
