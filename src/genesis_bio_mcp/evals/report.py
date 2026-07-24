"""Result types + Markdown/JSON rendering for the eval harness.

The report's **Capability gaps** section is the point of the whole exercise: it
lists the questions no current tool can answer (probes that failed) plus any
non-probe failures to investigate — i.e. the evidence for what to build next.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, field

# Two-sided 95% Student-t critical values by degrees of freedom (n-1). Hardcoded so the
# harness needs no scipy; df > 30 falls back to the normal approximation (1.96).
_T95 = {
    1: 12.706,
    2: 4.303,
    3: 3.182,
    4: 2.776,
    5: 2.571,
    6: 2.447,
    7: 2.365,
    8: 2.306,
    9: 2.262,
    10: 2.228,
    11: 2.201,
    12: 2.179,
    13: 2.160,
    14: 2.145,
    15: 2.131,
    16: 2.120,
    17: 2.110,
    18: 2.101,
    19: 2.093,
    20: 2.086,
    21: 2.080,
    22: 2.074,
    23: 2.069,
    24: 2.064,
    25: 2.060,
    26: 2.056,
    27: 2.052,
    28: 2.048,
    29: 2.045,
    30: 2.042,
}


def _mean_ci(
    values: list[float], *, bounds: tuple[float, float] | None = (0.0, 1.0)
) -> tuple[float, float | None, float | None]:
    """Mean and a t-based 95% CI over per-item values (the item is the sampling unit).

    Returns ``(mean, lo, hi)``; ``lo``/``hi`` are ``None`` when n < 2 (no spread to estimate).
    ``bounds`` clamps the interval to the quantity's natural range — ``(0, 1)`` for a rate,
    ``(-1, 1)`` for an uplift (difference of rates); pass ``None`` to leave it unclamped.
    """
    n = len(values)
    if n == 0:
        return 0.0, None, None
    mean = sum(values) / n
    if n < 2:
        return mean, None, None
    var = sum((v - mean) ** 2 for v in values) / (n - 1)
    half = _T95.get(n - 1, 1.96) * math.sqrt(var / n)
    lo, hi = mean - half, mean + half
    if bounds is not None:
        lo, hi = max(bounds[0], lo), min(bounds[1], hi)
    return mean, lo, hi


@dataclass
class ItemResult:
    """Outcome of one eval item across the deterministic (and optional agentic) layers."""

    id: str
    category: str
    probe: bool
    question: str

    skipped: bool = False
    skip_reason: str = ""

    # Deterministic layer (call declared tools directly, grade the output)
    det_passed: bool = False
    det_detail: str = ""
    det_latency_s: float = 0.0
    det_error: str = ""

    # Agentic layer (run the workflow agent end-to-end) — optional
    agentic_ran: bool = False
    agentic_passed: bool = False
    agentic_detail: str = ""
    tool_selection_ok: bool | None = None
    called_tools: list[str] = field(default_factory=list)

    # Prior-only baseline (same agent, NO tools) — optional; grades the tools-off answer
    baseline_ran: bool = False
    baseline_passed: bool = False
    baseline_detail: str = ""

    # Replicates (K): how many times the LLM layers ran, and how many of those passed. The
    # per-item rate is pass_count / replicates; the report aggregates rates across items.
    replicates: int = 1
    agentic_pass_count: int = 0
    baseline_pass_count: int = 0
    tool_selection_pass_count: int = 0


@dataclass
class EvalReport:
    generated_at: str
    results: list[ItemResult]
    agentic: bool
    baseline: bool = False
    replicates: int = 1

    # -- derived views -------------------------------------------------------
    @property
    def ran(self) -> list[ItemResult]:
        return [r for r in self.results if not r.skipped]

    @property
    def real(self) -> list[ItemResult]:
        return [r for r in self.ran if not r.probe]

    @property
    def probes(self) -> list[ItemResult]:
        return [r for r in self.ran if r.probe]

    @property
    def gaps(self) -> list[ItemResult]:
        """Probes the current tools could not answer — genuine capability gaps."""
        return [r for r in self.probes if not r.det_passed]

    @property
    def regressions(self) -> list[ItemResult]:
        """Non-probe failures — a miscalibrated grader or a real data-layer regression."""
        return [r for r in self.real if not r.det_passed]

    @property
    def shortcuts(self) -> list[ItemResult]:
        """Real items the prior-only (tools-off) baseline already answered — candidates to
        rewrite or cut, because they measure what the model knows, not what the tools add."""
        return [r for r in self.real if r.baseline_ran and r.baseline_passed]

    @property
    def tool_dependent(self) -> list[ItemResult]:
        """Real items only the tools made answerable: tools-on passed, tools-off failed."""
        return [
            r
            for r in self.real
            if r.agentic_ran and r.baseline_ran and r.agentic_passed and not r.baseline_passed
        ]

    def to_json_dict(self) -> dict:
        return {
            "generated_at": self.generated_at,
            "agentic": self.agentic,
            "replicates": self.replicates,
            "summary": self._summary(),
            "results": [
                {
                    "id": r.id,
                    "category": r.category,
                    "probe": r.probe,
                    "skipped": r.skipped,
                    "skip_reason": r.skip_reason,
                    "deterministic_passed": r.det_passed,
                    "deterministic_detail": r.det_detail,
                    "deterministic_latency_s": round(r.det_latency_s, 2),
                    "deterministic_error": r.det_error,
                    "agentic_ran": r.agentic_ran,
                    "agentic_passed": r.agentic_passed,
                    "tool_selection_ok": r.tool_selection_ok,
                    "called_tools": r.called_tools,
                    "baseline_ran": r.baseline_ran,
                    "baseline_passed": r.baseline_passed,
                    "replicates": r.replicates,
                    "agentic_pass_count": r.agentic_pass_count,
                    "baseline_pass_count": r.baseline_pass_count,
                    "tool_selection_pass_count": r.tool_selection_pass_count,
                }
                for r in self.results
            ],
        }

    def _summary(self) -> dict:
        real, probes = self.real, self.probes
        det_pass = sum(r.det_passed for r in real)
        s = {
            "deterministic_passed": det_pass,
            "deterministic_total": len(real),
            "probes_answerable": sum(r.det_passed for r in probes),
            "probes_total": len(probes),
            "capability_gaps": len(self.gaps),
            "skipped": sum(1 for r in self.results if r.skipped),
        }
        if self.agentic:
            agrun = [r for r in real if r.agentic_ran]
            sel = [r for r in agrun if r.tool_selection_ok is not None]
            s["agentic_answers_passed"] = sum(r.agentic_passed for r in agrun)
            s["agentic_total"] = len(agrun)
            s["tool_selection_ok"] = sum(bool(r.tool_selection_ok) for r in sel)
            s["tool_selection_total"] = len(sel)
        if self.baseline:
            both = [r for r in real if r.agentic_ran and r.baseline_ran]
            s["uplift_total"] = len(both)
            s["tools_on_passed"] = sum(r.agentic_passed for r in both)
            s["tools_off_passed"] = sum(r.baseline_passed for r in both)
            s["uplift_net"] = s["tools_on_passed"] - s["tools_off_passed"]
            s["tool_dependent"] = len(self.tool_dependent)
            s["shortcuts"] = len(self.shortcuts)
        # Replicate-aware confidence intervals over per-item rates (item = sampling unit).
        if self.replicates > 1 and self.agentic:

            def _ci(vals: list[float], bounds: tuple[float, float] = (0.0, 1.0)) -> dict:
                m, lo, hi = _mean_ci(vals, bounds=bounds)
                r3 = lambda x: None if x is None else round(x, 3)  # noqa: E731
                return {"mean": round(m, 3), "lo": r3(lo), "hi": r3(hi)}

            ag = [r for r in real if r.agentic_ran]
            s["replicates"] = self.replicates
            s["agentic_rate_ci"] = _ci([r.agentic_pass_count / r.replicates for r in ag])
            sel = [r for r in ag if r.tool_selection_ok is not None]
            if sel:
                s["tool_selection_rate_ci"] = _ci(
                    [r.tool_selection_pass_count / r.replicates for r in sel]
                )
            if self.baseline:
                both = [r for r in ag if r.baseline_ran]
                s["baseline_rate_ci"] = _ci([r.baseline_pass_count / r.replicates for r in both])
                s["uplift_ci"] = _ci(
                    [(r.agentic_pass_count - r.baseline_pass_count) / r.replicates for r in both],
                    bounds=(-1.0, 1.0),
                )
        return s

    # -- markdown ------------------------------------------------------------
    def to_markdown(self) -> str:
        s = self._summary()
        pct = (
            f"{100 * s['deterministic_passed'] / s['deterministic_total']:.0f}%"
            if s["deterministic_total"]
            else "n/a"
        )
        lines = [
            "# genesis-bio-mcp — evaluation report",
            "",
            f"_Generated {self.generated_at}_",
            "",
            "## Summary",
            f"- **Deterministic tasks:** {s['deterministic_passed']}/{s['deterministic_total']} passed ({pct})",
            f"- **Capability probes:** {s['probes_answerable']}/{s['probes_total']} answerable by current "
            f"tools → **{s['capability_gaps']} gaps**",
            f"- **Skipped:** {s['skipped']}",
        ]
        if self.agentic:
            sel_pct = (
                f"{100 * s['tool_selection_ok'] / s['tool_selection_total']:.0f}%"
                if s["tool_selection_total"]
                else "n/a"
            )
            lines.append(
                f"- **Agentic:** answers {s['agentic_answers_passed']}/{s['agentic_total']}; "
                f"tool-selection {s['tool_selection_ok']}/{s['tool_selection_total']} ({sel_pct})"
            )
        if self.baseline:
            lines.append(
                f"- **Uplift (do the tools help?):** tools-on {s['tools_on_passed']}/{s['uplift_total']} "
                f"vs tools-off {s['tools_off_passed']}/{s['uplift_total']} → **net {s['uplift_net']:+d}** "
                f"— {s['tool_dependent']} only answerable with tools, {s['shortcuts']} shortcuts pass without them"
            )
        if self.replicates > 1 and self.agentic:

            def _fmt(c: dict, pct: bool = True) -> str:
                if pct:
                    base = f"{100 * c['mean']:.0f}%"
                    return (
                        base
                        if c["lo"] is None
                        else f"{base} [{100 * c['lo']:.0f}, {100 * c['hi']:.0f}]"
                    )
                base = f"{c['mean']:+.2f}"
                return base if c["lo"] is None else f"{base} [{c['lo']:+.2f}, {c['hi']:+.2f}]"

            lines.append(
                f"- **Replicates:** K={self.replicates}/item; 95% CI is t-based across items "
                "(the item is the sampling unit — one run is an anecdote)"
            )
            lines.append(f"  - Agentic answer rate: {_fmt(s['agentic_rate_ci'])}")
            if "tool_selection_rate_ci" in s:
                lines.append(f"  - Tool-selection rate: {_fmt(s['tool_selection_rate_ci'])}")
            if "baseline_rate_ci" in s:
                lines.append(f"  - Tools-off baseline rate: {_fmt(s['baseline_rate_ci'])}")
            if "uplift_ci" in s:
                lines.append(f"  - Uplift (Δ pass-rate): {_fmt(s['uplift_ci'], pct=False)}")

        # Coverage by category (real items only)
        by_cat: dict[str, list[ItemResult]] = defaultdict(list)
        for r in self.real:
            by_cat[r.category].append(r)
        lines += ["", "## Coverage by category", "| Category | Passed | Total |", "|---|---:|---:|"]
        for cat in sorted(by_cat):
            rs = by_cat[cat]
            lines.append(f"| {cat} | {sum(x.det_passed for x in rs)} | {len(rs)} |")

        # Capability gaps — the headline output
        lines += ["", "## Capability gaps — what to build next"]
        if self.gaps:
            for r in self.gaps:
                lines.append(f"- **{r.id}** ({r.category}): {r.question}")
        else:
            lines.append("_No unmet probes — every capability probe was answerable._")
        if self.regressions:
            lines += [
                "",
                "### Non-probe failures (investigate — grader miscalibration or a real regression)",
            ]
            for r in self.regressions:
                why = r.det_error or r.det_detail
                lines.append(f"- **{r.id}** ({r.category}): {why}")

        # Uplift — the whole point of the baseline: what do the tools actually add?
        if self.baseline:
            td, sc = self.tool_dependent, self.shortcuts
            lines += ["", "## Tool-value analysis (uplift)"]
            lines.append(
                f"- **Only answerable with tools ({len(td)}):** "
                + (
                    ", ".join(r.id for r in td)
                    or "_none — the tools added nothing the model didn't already know_"
                )
            )
            lines.append(
                f"- **Shortcuts — pass tools-off, rewrite or cut ({len(sc)}):** "
                + (", ".join(r.id for r in sc) or "_none_")
            )

        # Per-item results
        lines += ["", "## Results"]
        if self.agentic:
            hdr = "| Item | Category | Deterministic | Latency (s) | Agentic | Tools ✓ |"
            sep = "|---|---|---|---:|---|---|"
            if self.baseline:
                hdr += " Tools-off |"
                sep += "---|"
            lines += [hdr, sep]
        else:
            lines += [
                "| Item | Category | Deterministic | Latency (s) |",
                "|---|---|---|---:|",
            ]
        for r in self.results:
            if r.skipped:
                det = f"skipped ({r.skip_reason})"
            elif r.probe:
                det = "GAP" if not r.det_passed else "✓ (answerable!)"
            else:
                det = "✓ pass" if r.det_passed else f"✗ {r.det_detail}"[:60]
            row = f"| {r.id} | {r.category} | {det} | {r.det_latency_s:.2f} |"
            if self.agentic:
                ag = "—" if not r.agentic_ran else ("✓" if r.agentic_passed else "✗")
                tools = (
                    "—" if r.tool_selection_ok is None else ("✓" if r.tool_selection_ok else "✗")
                )
                row += f" {ag} | {tools} |"
                if self.baseline:
                    off = "—" if not r.baseline_ran else ("✓" if r.baseline_passed else "✗")
                    row += f" {off} |"
            lines.append(row)

        return "\n".join(lines)
