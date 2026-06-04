"""Pure grading functions for the eval harness.

No I/O, no LLM — each grader maps (output_text, spec) → GradeResult, so they are
fully unit-testable and run in CI. A grader spec is a small dict with a ``type``:

- ``contains``  — case-insensitive substring: ``{"type": "contains", "value": "P04626"}``
- ``all_of``    — every substring present: ``{"type": "all_of", "value": ["a", "b"]}``
- ``regex``     — pattern search (ignorecase by default): ``{"type": "regex", "value": "p.?Chembl\\s*9"}``
- ``numeric``   — a number near a label satisfies a bound:
                  ``{"type": "numeric", "near": "pChEMBL", "min": 9.0}`` or
                  ``{"type": "numeric", "value": 180.16, "tolerance": 1.0}``
"""

from __future__ import annotations

import re
from dataclasses import dataclass

_NUM_RE = re.compile(r"[-+]?\d+(?:\.\d+)?")


@dataclass
class GradeResult:
    """Outcome of applying one grader to one output string."""

    passed: bool
    detail: str


def grade(output: str | None, spec: dict) -> GradeResult:
    """Apply a grader spec to a tool/agent output string. Never raises."""
    text = output or ""
    gtype = spec.get("type")
    dispatch = {
        "contains": _contains,
        "all_of": _all_of,
        "regex": _regex,
        "numeric": _numeric,
    }
    fn = dispatch.get(gtype)
    if fn is None:
        return GradeResult(False, f"unknown grader type: {gtype!r}")
    try:
        return fn(text, spec)
    except Exception as exc:  # a malformed spec should fail the item, not crash the run
        return GradeResult(False, f"grader error: {exc}")


def _contains(text: str, spec: dict) -> GradeResult:
    needle = str(spec["value"])
    ok = needle.lower() in text.lower()
    return GradeResult(ok, f"{'found' if ok else 'missing'} {needle!r}")


def _all_of(text: str, spec: dict) -> GradeResult:
    needles = [str(v) for v in spec["value"]]
    low = text.lower()
    missing = [n for n in needles if n.lower() not in low]
    return GradeResult(not missing, "all present" if not missing else f"missing {missing}")


def _regex(text: str, spec: dict) -> GradeResult:
    pattern = str(spec["value"])
    flags = 0 if spec.get("ignorecase") is False else re.IGNORECASE
    ok = re.search(pattern, text, flags) is not None
    return GradeResult(ok, f"{'matched' if ok else 'no match'} /{pattern}/")


def _numeric(text: str, spec: dict) -> GradeResult:
    """A number (optionally near a label) must satisfy value±tol or min/max bounds."""
    window = text
    near = spec.get("near")
    if near:
        idx = text.lower().find(str(near).lower())
        if idx == -1:
            return GradeResult(False, f"label {near!r} not found")
        window = text[idx : idx + 80]
    nums = [float(m) for m in _NUM_RE.findall(window)]
    if not nums:
        return GradeResult(False, "no number found")

    if "value" in spec:
        expected = float(spec["value"])
        tol = float(spec.get("tolerance", 0.0))
        hit = next((n for n in nums if abs(n - expected) <= tol), None)
        if hit is None:
            return GradeResult(False, f"no number within {tol} of {expected} (saw {nums[:8]})")
        return GradeResult(True, f"{hit} ≈ {expected}")

    lo = float(spec["min"]) if "min" in spec else float("-inf")
    hi = float(spec["max"]) if "max" in spec else float("inf")
    hit = next((n for n in nums if lo <= n <= hi), None)
    if hit is None:
        return GradeResult(False, f"no number in [{lo}, {hi}] (saw {nums[:8]})")
    return GradeResult(True, f"{hit} in [{lo}, {hi}]")
