"""CLI entrypoint for the eval harness.

    uv run python -m genesis_bio_mcp.evals [--agentic] [--out report.md] [--filter CATEGORY]

Always runs the deterministic (keyless) layer; ``--agentic`` adds the workflow-agent
layer when ANTHROPIC_API_KEY is set. Prints the Markdown report; ``--out`` also writes
it plus a sibling ``.json``.
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
    report = asyncio.run(
        run_eval(
            items, agentic=args.agentic, max_iterations=args.max_iterations, timestamp=timestamp
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
