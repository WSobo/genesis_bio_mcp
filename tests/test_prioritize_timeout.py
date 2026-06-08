"""Regression tests for the compare_targets / prioritize_target wedging fix.

A slow upstream source used to hang a whole gene assessment indefinitely (no
per-sub-query timeout), and a 5-gene compare_targets fan-out could leave the
server wedged: the orphaned work kept holding the small per-source semaphores so
every later call starved. These tests lock in the fix:

- ``_safe`` bounds each sub-query with ``asyncio.wait_for`` and degrades to a
  "timed out" data gap instead of hanging;
- a timed-out sub-query *releases* the client semaphore it held (the actual
  anti-wedge guarantee);
- ``prioritize_target`` returns promptly with the slow source marked as an error;
- the GWAS study fan-out and Open Targets disease-variant loop are capped.
"""

from __future__ import annotations

import asyncio
from unittest.mock import AsyncMock, MagicMock

import httpx

import genesis_bio_mcp.clients.open_targets as ot_mod
import genesis_bio_mcp.tools.target_prioritization as tp
from genesis_bio_mcp.clients.gwas import GwasClient
from genesis_bio_mcp.clients.open_targets import OpenTargetsClient
from genesis_bio_mcp.config.settings import settings
from genesis_bio_mcp.models import GeneResolution
from genesis_bio_mcp.tools.target_prioritization import _safe, _safe_timed, prioritize_target


async def test_safe_returns_result_on_success() -> None:
    async def ok():
        return 42

    result, err = await _safe(ok(), timeout=5.0)
    assert result == 42
    assert err is None


async def test_safe_preserves_exception_behavior() -> None:
    async def boom():
        raise ValueError("nope")

    result, err = await _safe(boom(), fallback="fb", timeout=5.0)
    assert result == "fb"
    assert "nope" in err


async def test_safe_times_out_slow_coro_quickly() -> None:
    async def slow():
        await asyncio.sleep(5)
        return "never"

    t0 = asyncio.get_running_loop().time()
    result, err = await _safe(slow(), fallback=None, timeout=0.1)
    elapsed = asyncio.get_running_loop().time() - t0

    assert result is None
    assert err is not None and "timed out" in err
    assert elapsed < 1.0  # bounded, did not wait the full 5s


async def test_safe_timed_passes_timeout_through() -> None:
    async def slow():
        await asyncio.sleep(5)

    _, err, secs = await _safe_timed("x", slow(), timeout=0.1)
    assert err is not None and "timed out" in err
    assert secs < 1.0


async def test_timed_out_subquery_releases_semaphore() -> None:
    """The anti-wedge guarantee: cancelling a hung sub-query frees the semaphore
    it held, so subsequent work can't starve."""
    sem = asyncio.Semaphore(1)

    async def hog():
        async with sem:
            await asyncio.sleep(5)

    result, err = await _safe(hog(), timeout=0.1)
    assert result is None
    assert "timed out" in err

    # The cancelled coroutine must have released the semaphore.
    assert sem.locked() is False
    # And it is immediately re-acquirable — no leak, no wedge.
    await asyncio.wait_for(sem.acquire(), timeout=0.5)
    sem.release()


async def test_prioritize_target_marks_slow_source_as_error(monkeypatch) -> None:
    """End-to-end: one hung source becomes a timed-out data gap; the report still
    returns promptly instead of hanging the whole assessment."""
    monkeypatch.setattr(settings, "prioritize_source_budget_secs", 0.1)
    monkeypatch.setattr(
        tp,
        "resolve_gene",
        AsyncMock(return_value=GeneResolution(hgnc_symbol="TP53", source="test")),
    )

    async def _hang(*_args, **_kwargs):
        await asyncio.sleep(5)
        return None

    uniprot = AsyncMock()
    uniprot.get_protein = AsyncMock(return_value=None)
    open_targets = AsyncMock()
    open_targets.get_association = AsyncMock(return_value=None)
    depmap = AsyncMock()
    depmap.get_essentiality = AsyncMock(return_value=None)
    gwas = AsyncMock()
    gwas.get_evidence = AsyncMock(return_value=None)
    pubchem = AsyncMock()
    pubchem.get_compounds = AsyncMock(return_value=None)
    chembl = AsyncMock()
    chembl.get_compounds = _hang  # the slow source

    t0 = asyncio.get_running_loop().time()
    report = await prioritize_target(
        "TP53",
        "cancer",
        uniprot=uniprot,
        open_targets=open_targets,
        depmap=depmap,
        gwas=gwas,
        pubchem=pubchem,
        chembl=chembl,
    )
    elapsed = asyncio.get_running_loop().time() - t0

    assert elapsed < 2.0  # did not hang on the 5s source
    assert "chembl" in report.data_gaps
    assert "timed out" in report.errors.get("chembl", "")


async def test_gwas_study_fanout_is_capped(monkeypatch) -> None:
    monkeypatch.setattr(settings, "gwas_study_fanout_cap", 2)
    resp = MagicMock()
    resp.json.return_value = {}
    resp.raise_for_status = MagicMock()

    async with httpx.AsyncClient() as http:
        gwas = GwasClient(http)
        gwas._client.get = AsyncMock(return_value=resp)
        assocs = [{"_links": {"study": {"href": f"http://gwas/{i}/study"}}} for i in range(5)]
        await gwas._resolve_study_data(assocs)

    # 5 distinct study links, but only `cap` sub-fetches are issued.
    assert gwas._client.get.call_count == 2


async def test_open_targets_variant_fanout_is_capped(monkeypatch) -> None:
    monkeypatch.setattr(settings, "open_targets_disease_variant_cap", 3)
    monkeypatch.setattr(
        ot_mod, "_normalize_indication_variants", lambda name: [f"v{i}" for i in range(10)]
    )

    async with httpx.AsyncClient() as http:
        ot = OpenTargetsClient(http)
        graphql = AsyncMock(return_value=None)
        monkeypatch.setattr(ot, "_graphql", graphql)
        result = await ot._resolve_disease("some disease")

    # 10 variants available, but disease resolution stops after `cap` GraphQL calls.
    assert graphql.call_count == 3
    assert result == (None, None)
