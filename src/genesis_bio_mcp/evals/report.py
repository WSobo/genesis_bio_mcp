"""Result types + Markdown/JSON rendering for the eval harness.

The report's **Capability gaps** section is the point of the whole exercise: it
lists the questions no current tool can answer (probes that failed) plus any
non-probe failures to investigate — i.e. the evidence for what to build next.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field


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


@dataclass
class EvalReport:
    generated_at: str
    results: list[ItemResult]
    agentic: bool

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

    def to_json_dict(self) -> dict:
        return {
            "generated_at": self.generated_at,
            "agentic": self.agentic,
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

        # Per-item results
        lines += ["", "## Results"]
        if self.agentic:
            lines += [
                "| Item | Category | Deterministic | Latency (s) | Agentic | Tools ✓ |",
                "|---|---|---|---:|---|---|",
            ]
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
                row = row[:-1] + f" {ag} | {tools} |"
            lines.append(row)

        return "\n".join(lines)
