"""Pure grading functions for the eval harness.

No I/O, no LLM — each grader maps (output_text, spec) → GradeResult, so they are
fully unit-testable and run in CI. A grader spec is a small dict with a ``type``:

- ``contains``  — case-insensitive substring: ``{"type": "contains", "value": "P04626"}``
- ``all_of``    — every substring present: ``{"type": "all_of", "value": ["a", "b"]}``
- ``regex``     — pattern search (ignorecase by default): ``{"type": "regex", "value": "p.?Chembl\\s*9"}``
- ``numeric``   — a number near a label satisfies a bound:
                  ``{"type": "numeric", "near": "pChEMBL", "min": 9.0}`` or
                  ``{"type": "numeric", "value": 180.16, "tolerance": 1.0}``
- ``set``       — recall (and optionally precision) over a canonical set. The right grader for
                  list answers (off-target kinases, similar compounds, pathway/ID sets) —
                  generalizes ``all_of`` (which is ``recall_at_k`` = 1.0) and is far stronger than
                  ``contains``. Three layers, weakest to strongest:
                    * recall-only (substring):
                      ``{"type": "set", "canonical": ["ABL1", "LYN", "FYN"], "recall_at_k": 0.66}``
                    * recall + **precision** — add an ``extractor`` regex whose whole matches (capturing
                      groups are ignored) delimit the *predicted* set from the free markdown — the answer
                      the run committed to — so over-listing every candidate can no longer game recall:
                      ``{"type": "set", "canonical": ["R-HSA-187706", ...], "extractor": "R-HSA-\\d+",
                         "recall_at_k": 0.6, "precision_at_k": 0.5}``
                    * ``forbidden`` — any listed distractor present (substring) fails the item outright;
                      combine with the above to punish the specific wrong answer a naive path emits.
                  This realises the mcp-eval-authoring skill's ``marker_gene_precision_recall`` P@K/R@K
                  semantics under the existing ``set`` type. ``precision_at_k`` REQUIRES an ``extractor``
                  (precision is undefined without a delimited predicted set).
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
        "set": _set,
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


def _set(text: str, spec: dict) -> GradeResult:
    """Recall — and optionally precision — over a canonical set (see module docstring).

    Without an ``extractor`` this is recall-only, case-insensitive **substring** membership: the
    fraction of ``canonical`` present in the output must meet ``recall_at_k`` (default 1.0). This
    is the historical behaviour and is preserved byte-for-byte.

    With an ``extractor`` regex, its matches form the *predicted* set (the answer the run actually
    committed to), membership becomes exact case-insensitive token equality, and ``precision_at_k``
    becomes checkable — recall alone can always be gamed by listing every candidate, precision
    cannot. ``forbidden`` fails outright on any listed distractor, with or without an extractor.
    """
    canonical = [str(v) for v in spec["canonical"]]
    if not canonical:
        return GradeResult(False, "empty canonical set")
    low = text.lower()

    forbidden = [str(v) for v in spec.get("forbidden", [])]
    recall_k = float(spec.get("recall_at_k", 1.0))
    precision_k = spec.get("precision_at_k")
    extractor = spec.get("extractor")

    if extractor is None:
        if precision_k is not None:
            return GradeResult(
                False, "precision_at_k requires an 'extractor' regex to delimit the predicted set"
            )
        # No predicted set to check against, so 'forbidden' is a coarse whole-text substring guard.
        hit_forbidden = [f for f in forbidden if f.lower() in low]
        if hit_forbidden:
            return GradeResult(False, f"forbidden token(s) present: {hit_forbidden[:8]}")
        present = [c for c in canonical if c.lower() in low]
        recall = len(present) / len(canonical)
        missing = [c for c in canonical if c not in present]
        detail = f"recall {len(present)}/{len(canonical)} = {recall:.2f} (need ≥ {recall_k})"
        if missing:
            detail += f"; missing {missing[:8]}"
        return GradeResult(recall >= recall_k, detail)

    # Extractor path: the regex matches ARE the predicted set, so precision is well-defined.
    # ``group(0)`` (via finditer) so a capturing group in the pattern can't corrupt the token set;
    # both sides of recall are set-based, so a duplicate ``canonical`` entry can't skew the ratio.
    pred_lower = {m.group(0).lower() for m in re.finditer(str(extractor), text, re.IGNORECASE)}
    canon_lower = {c.lower() for c in canonical}
    # Here 'forbidden' is checked against the committed answer tokens (exact), matching precision.
    hit_forbidden = [f for f in forbidden if f.lower() in pred_lower]
    if hit_forbidden:
        return GradeResult(False, f"forbidden token(s) present: {hit_forbidden[:8]}")
    hits = canon_lower & pred_lower
    recall = len(hits) / len(canon_lower)
    precision = len(hits) / len(pred_lower) if pred_lower else 0.0

    ok = recall >= recall_k
    detail = f"recall {len(hits)}/{len(canon_lower)} = {recall:.2f} (≥ {recall_k})"
    if precision_k is not None:
        precision_k = float(precision_k)
        ok = ok and precision >= precision_k
        detail += f", precision {len(hits)}/{len(pred_lower)} = {precision:.2f} (≥ {precision_k})"
    else:
        detail += f", precision {precision:.2f} ({len(pred_lower)} predicted)"
    missing = [c for c in canonical if c.lower() not in pred_lower]
    if missing:
        detail += f"; missing {missing[:8]}"
    return GradeResult(ok, detail)


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
