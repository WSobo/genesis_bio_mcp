"""IEDB NextGen Tools API client — MHC-I and MHC-II binding prediction.

The NextGen Tools service (``https://api-nextgen-tools.iedb.org/api/v1/``)
is the modern HTTPS replacement for ``http://tools-cluster-interface.iedb.org``.
It uses an async job pattern: POST a pipeline spec, poll the returned
``results_uri`` until ``status == "done"``, then read the results table.

Methods available (as of 2026-04-16, confirmed via the ``/mhci`` self-docs
endpoint):
  * ``netmhcpan_el``  — IEDB-recommended epitope predictor (NetMHCpan 4.1 EL)
  * ``netmhcpan_ba``  — IEDB-recommended binding-affinity predictor
  * ``netmhciipan_el`` / ``netmhciipan_ba`` — class-II equivalents
  * ``mhcflurry``, ``smm``, ``smmpmbec``, ``ann``, ``consensus``, ...

Only the EL-mode methods are exposed as defaults; the others are reachable
via the optional ``method`` parameter.
"""

from __future__ import annotations

import asyncio
import logging
from typing import Any, Literal

import httpx

from genesis_bio_mcp.models import MHCBindingHit, MHCBindingResults

logger = logging.getLogger(__name__)

_BASE_URL = "https://api-nextgen-tools.iedb.org/api/v1"
_PIPELINE_URL = f"{_BASE_URL}/pipeline"
_SEMAPHORE = asyncio.Semaphore(2)  # shared academic compute — be polite

# IEDB-convention binder thresholds.
# Source: IEDB Analysis Tools documentation (tools.iedb.org/mhci/help/).
# Tuning note: loosening to 1.0 (strong) / 5.0 (weak) is common when
# screening large peptide pools for secondary analysis; the canonical
# 0.5 / 2.0 pair is what IEDB reports in peer-reviewed benchmarks.
_STRONG_BINDER_PCTILE: float = 0.5
_WEAK_BINDER_PCTILE: float = 2.0

# Default HLA panels. Coverage figures are IEDB/literature consensus,
# not computed from our own population files.
#
# Class I: 5-allele IEDB reference set, ~85% global phenotype coverage.
# Source: Weiskopf et al., J Immunol 2013 (IEDB class-I reference).
_DEFAULT_HLA_CLASS_I: tuple[str, ...] = (
    "HLA-A*02:01",
    "HLA-A*03:01",
    "HLA-A*24:02",
    "HLA-B*07:02",
    "HLA-B*35:01",
)
# Class II: 5-allele DRB1 panel, ~60% global coverage.
# Source: Greenbaum et al., Immunogenetics 2011.
_DEFAULT_HLA_CLASS_II: tuple[str, ...] = (
    "HLA-DRB1*01:01",
    "HLA-DRB1*03:01",
    "HLA-DRB1*04:01",
    "HLA-DRB1*07:01",
    "HLA-DRB1*15:01",
)

# Canonical peptide lengths for each MHC class per IEDB best practice.
_DEFAULT_LENGTHS_CLASS_I: tuple[int, int] = (9, 10)
_DEFAULT_LENGTHS_CLASS_II: tuple[int, int] = (15, 15)

# Default predictor methods (EL = Elution Ligand mode, IEDB-recommended
# since 2023.09 per the /mhci self-docs endpoint).
_DEFAULT_METHOD_CLASS_I: str = "netmhcpan_el"
_DEFAULT_METHOD_CLASS_II: str = "netmhciipan_el"

# Polling parameters. The NextGen service typically completes a small job
# in 2–10s; we cap at 60s total wall time to bound MCP tool latency.
#
# Tuning note: raising _POLL_TIMEOUT_SEC lets larger panels (e.g. 20+ HLAs
# × 100 peptides) complete, at the cost of blocking the tool for up to
# that many seconds. Prefer narrowing the input (smaller panel, shorter
# peptide range) over raising the timeout.
_POLL_INTERVAL_SEC: float = 2.0
_POLL_TIMEOUT_SEC: float = 60.0

# Request-size cap. IEDB's shared service throttles heavy jobs; exceeding
# this is almost always a user error (feed a whole protein + 20 alleles).
_MAX_REQUEST_SIZE: int = 2000  # approx peptides × alleles


class IEDBToolsClient:
    """Client for IEDB NextGen Tools MHC-I / MHC-II binding prediction."""

    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client

    async def predict_mhc_binding(
        self,
        sequence: str,
        alleles: list[str] | None = None,
        mhc_class: Literal["I", "II"] = "I",
        peptide_lengths: list[int] | None = None,
        method: str | None = None,
    ) -> MHCBindingResults | None:
        """Submit a binding prediction job and poll for completion.

        Args:
            sequence: Peptide, FASTA, or full protein. If not already FASTA,
                wrapped as ``">query\\n{sequence}"``. The IEDB service
                handles peptide windowing internally based on ``peptide_lengths``.
            alleles: HLA allele list. Defaults to the 5-allele class-I or
                class-II reference panel (see module constants).
            mhc_class: ``"I"`` (default) or ``"II"``.
            peptide_lengths: 2-list ``[min, max]`` of peptide window sizes.
                Defaults to ``[9, 10]`` for class I, ``[15, 15]`` for class II.
            method: Predictor method. Defaults to ``netmhcpan_el`` (class I)
                or ``netmhciipan_el`` (class II).

        Returns:
            :class:`MHCBindingResults` on success, ``None`` on network error
            or API rejection. Times out with a partial result and a note if
            the job doesn't finish within :data:`_POLL_TIMEOUT_SEC`.

        Raises:
            ValueError: if the estimated peptide-count × allele-count
                product exceeds :data:`_MAX_REQUEST_SIZE`.
        """
        if mhc_class not in ("I", "II"):
            raise ValueError(f"mhc_class must be 'I' or 'II', got {mhc_class!r}")

        panel = (
            list(alleles)
            if alleles
            else list(_DEFAULT_HLA_CLASS_I if mhc_class == "I" else _DEFAULT_HLA_CLASS_II)
        )
        if not panel:
            raise ValueError("alleles list is empty")

        default_range = _DEFAULT_LENGTHS_CLASS_I if mhc_class == "I" else _DEFAULT_LENGTHS_CLASS_II
        if peptide_lengths:
            if len(peptide_lengths) == 1:
                length_range = [peptide_lengths[0], peptide_lengths[0]]
            else:
                length_range = [min(peptide_lengths), max(peptide_lengths)]
        else:
            length_range = list(default_range)

        method = method or (
            _DEFAULT_METHOD_CLASS_I if mhc_class == "I" else _DEFAULT_METHOD_CLASS_II
        )

        fasta_text = _ensure_fasta(sequence)
        peptide_count_est = _estimate_peptide_count(fasta_text, length_range)
        request_size = peptide_count_est * len(panel)
        if request_size > _MAX_REQUEST_SIZE:
            raise ValueError(
                f"Request size ({peptide_count_est} peptide windows × "
                f"{len(panel)} alleles = {request_size}) exceeds limit of "
                f"{_MAX_REQUEST_SIZE}. Narrow the input sequence, reduce the "
                "HLA panel, or tighten peptide_lengths."
            )

        payload = _build_payload(mhc_class, fasta_text, panel, length_range, method)

        async with _SEMAPHORE:
            submit = await self._submit(payload)
            if submit is None:
                return None
            results_uri = submit.get("results_uri")
            if not results_uri:
                logger.warning("IEDB submit response missing results_uri: %s", submit)
                return None
            poll_result, timed_out = await self._poll(results_uri)

        if poll_result is None:
            return None

        hits = _parse_results(poll_result, mhc_class=mhc_class, method=method)
        strong = sum(
            1
            for h in hits
            if h.percentile_rank is not None and h.percentile_rank < _STRONG_BINDER_PCTILE
        )
        weak = sum(
            1
            for h in hits
            if h.percentile_rank is not None
            and _STRONG_BINDER_PCTILE <= h.percentile_rank < _WEAK_BINDER_PCTILE
        )
        notes: list[str] = []
        if timed_out:
            notes.append(
                f"IEDB job did not finish within {_POLL_TIMEOUT_SEC:.0f}s; "
                "results may be incomplete — retry with a narrower panel or shorter sequence."
            )
        return MHCBindingResults(
            input_sequence=sequence.strip(),
            mhc_class=mhc_class,
            method=method,
            alleles_tested=panel,
            peptide_length_range=(length_range[0], length_range[1]),
            hits=hits,
            strong_binder_count=strong,
            weak_binder_count=weak,
            notes=notes,
        )

    async def _submit(self, payload: dict) -> dict | None:
        try:
            resp = await self._client.post(_PIPELINE_URL, json=payload, timeout=30.0)
            resp.raise_for_status()
            return resp.json()
        except Exception as exc:
            logger.warning("IEDB submit failed: %s", exc)
            return None

    async def _poll(self, uri: str) -> tuple[dict | None, bool]:
        """Poll the results URI until status='done' or timeout.

        Returns (result_payload, timed_out). If the service reports ``done``
        partway through, ``timed_out`` is False and the payload is the final
        result. If the loop exits because the deadline elapsed, we return
        the last response and ``timed_out=True``.
        """
        elapsed = 0.0
        last: dict | None = None
        while elapsed < _POLL_TIMEOUT_SEC:
            try:
                resp = await self._client.get(uri, timeout=20.0)
                resp.raise_for_status()
                last = resp.json()
            except Exception as exc:
                logger.warning("IEDB poll failed: %s", exc)
                return None, False
            if isinstance(last, dict) and last.get("status") == "done":
                return last, False
            await asyncio.sleep(_POLL_INTERVAL_SEC)
            elapsed += _POLL_INTERVAL_SEC
        return last, True


# ---------------------------------------------------------------------------
# Payload + parse helpers
# ---------------------------------------------------------------------------


def _ensure_fasta(seq: str) -> str:
    """Wrap a bare sequence in a FASTA header if it doesn't already have one."""
    s = seq.strip()
    if s.startswith(">"):
        return s
    return f">query\n{s}"


def _estimate_peptide_count(fasta: str, length_range: list[int]) -> int:
    """Estimate the number of windowed peptides the IEDB service will produce."""
    count = 0
    min_len = length_range[0]
    max_len = length_range[1]
    for line in fasta.splitlines():
        if line.startswith(">"):
            continue
        L = len(line.strip())
        if L == 0:
            continue
        if L < min_len:
            continue
        # Rough count: for each length l in [min_len, max_len], (L - l + 1) windows.
        for length in range(min_len, max_len + 1):
            n = L - length + 1
            if n > 0:
                count += n
    return max(count, 1)


def _build_payload(
    mhc_class: str,
    fasta_text: str,
    alleles: list[str],
    length_range: list[int],
    method: str,
) -> dict:
    """Build the NextGen /pipeline POST body for a single-stage MHC job."""
    return {
        "pipeline_id": "",
        "run_stage_range": [1, 1],
        "stages": [
            {
                "stage_number": 1,
                "tool_group": "mhci" if mhc_class == "I" else "mhcii",
                "input_sequence_text": fasta_text,
                "input_parameters": {
                    "alleles": ",".join(alleles),
                    "peptide_length_range": length_range,
                    "predictors": [{"type": "binding", "method": method}],
                },
            }
        ],
    }


def _parse_results(
    payload: dict, *, mhc_class: Literal["I", "II"], method: str
) -> list[MHCBindingHit]:
    """Extract MHCBindingHit rows from a completed NextGen result payload.

    The payload shape (verified live against the /mhci endpoint on
    2026-04-16) is::

        {
          "status": "done",
          "data": {
            "results": [
              {"type": "peptide_table",
               "table_columns": [{"name": "peptide", ...}, ...],
               "table_data": [[val, val, ...], ...]},
              {"type": "netmhcpan_allele_distance", ...},
              ...
            ]
          }
        }

    We look for the first ``peptide_table`` block and zip its rows against
    the column metadata. Column names differ per method — we look up
    ``peptide``, ``allele``, ``length``, ``median_percentile``, and the
    method-specific ``{method}_percentile`` / ``{method}_score`` /
    ``{method}_core`` columns.
    """
    results = (payload.get("data") or {}).get("results") or []
    peptide_table: dict | None = None
    for block in results:
        if isinstance(block, dict) and block.get("type") == "peptide_table":
            peptide_table = block
            break
    if peptide_table is None:
        return []
    columns = peptide_table.get("table_columns") or []
    rows = peptide_table.get("table_data") or []
    col_index: dict[str, int] = {
        c.get("name"): i for i, c in enumerate(columns) if isinstance(c, dict)
    }

    pct_key = f"{method}_percentile"
    score_key = f"{method}_score"
    core_key = f"{method}_core"

    hits: list[MHCBindingHit] = []
    for row in rows:
        if not isinstance(row, list):
            continue
        pct = _get(row, col_index, pct_key)
        if pct is None:
            pct = _get(row, col_index, "median_percentile")
        score = _get(row, col_index, score_key)
        core = _get(row, col_index, core_key)
        peptide = _get(row, col_index, "peptide")
        allele = _get(row, col_index, "allele")
        length = _get(row, col_index, "length")
        if peptide is None or allele is None:
            continue
        try:
            pct_f = float(pct) if pct is not None else None
        except (TypeError, ValueError):
            pct_f = None
        try:
            score_f = float(score) if score is not None else None
        except (TypeError, ValueError):
            score_f = None
        try:
            length_i = int(length) if length is not None else None
        except (TypeError, ValueError):
            length_i = None

        if pct_f is None:
            binder_class: Literal["strong", "weak", "non_binder"] = "non_binder"
        elif pct_f < _STRONG_BINDER_PCTILE:
            binder_class = "strong"
        elif pct_f < _WEAK_BINDER_PCTILE:
            binder_class = "weak"
        else:
            binder_class = "non_binder"

        hits.append(
            MHCBindingHit(
                peptide=str(peptide),
                allele=str(allele),
                peptide_length=length_i,
                percentile_rank=pct_f,
                score=score_f,
                core_peptide=str(core) if core else None,
                binder_class=binder_class,
            )
        )
    # Sort strongest binder (lowest percentile) first
    hits.sort(key=lambda h: h.percentile_rank if h.percentile_rank is not None else 999.0)
    return hits


def _get(row: list, col_index: dict[str, int], name: str) -> Any:
    idx = col_index.get(name)
    if idx is None or idx >= len(row):
        return None
    return row[idx]
