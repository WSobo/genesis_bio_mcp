"""Unit tests for the v0.6.0 corpus store scaffold (M0).

These cover the DB-optional behavior and the corpus_describe tool WITHOUT a live database
(the corpus tier is optional by design). A real-Postgres integration test lives separately and
is skipped unless GENESIS_CORPUS_DSN is set.
"""

from __future__ import annotations

import json
import os
from datetime import UTC, datetime
from types import SimpleNamespace

import pytest

from genesis_bio_mcp import server
from genesis_bio_mcp.corpus import create_corpus_pool, fetch_manifest
from genesis_bio_mcp.corpus import db as corpus_db
from genesis_bio_mcp.models import CorpusManifest
from genesis_bio_mcp.server import CorpusDescribeInput, corpus_describe


def test_schema_sql_is_bundled_and_nonempty():
    sql = corpus_db._schema_sql()
    assert "corpus_manifest" in sql
    assert "CREATE EXTENSION IF NOT EXISTS vector" in sql
    # All four tables are declared.
    for table in ("targets", "compounds", "activities", "corpus_manifest"):
        assert table in sql


@pytest.mark.asyncio
async def test_create_corpus_pool_returns_none_without_dsn(monkeypatch):
    # No DSN configured → optional corpus tier yields None, never raises, never touches a DB.
    monkeypatch.setattr(corpus_db.settings, "corpus_dsn", None)
    assert await corpus_db.create_corpus_pool() is None


@pytest.mark.asyncio
async def test_corpus_describe_unconfigured_returns_upstream_unavailable():
    server.mcp.state = SimpleNamespace(corpus_pool=None)
    out = await corpus_describe(CorpusDescribeInput(response_format="json"))
    env = json.loads(out)
    assert env["error"]["status"] == "UpstreamUnavailable"
    assert "not configured" in env["error"]["message"]


@pytest.mark.asyncio
async def test_corpus_describe_unbuilt_returns_not_found(monkeypatch):
    # Pool present but the manifest row is absent → corpus configured but not yet built.
    server.mcp.state = SimpleNamespace(corpus_pool=object())

    async def _no_manifest(_pool):
        return None

    monkeypatch.setattr(server, "fetch_manifest", _no_manifest)
    out = await corpus_describe(CorpusDescribeInput(response_format="json"))
    env = json.loads(out)
    assert env["error"]["status"] == "NotFound"
    assert "not been built" in env["error"]["message"]


@pytest.mark.asyncio
async def test_corpus_describe_returns_manifest(monkeypatch):
    server.mcp.state = SimpleNamespace(corpus_pool=object())

    async def _manifest(_pool):
        return {
            "target_family": "human kinome",
            "built_at": datetime(2026, 5, 30, tzinfo=UTC),
            "chembl_release": "ChEMBL_34",
            "uniprot_snapshot": "2026_02",
            "protein_embedding_model": "ESM-2 650M (mean-pooled)",
            "chem_embedding_model": None,
            "target_count": 512,
            "compound_count": 41234,
            "activity_count": 87654,
        }

    monkeypatch.setattr(server, "fetch_manifest", _manifest)
    out = await corpus_describe(CorpusDescribeInput(response_format="json"))
    env = json.loads(out)
    assert env["provenance"]["source"] == "genesis-bio-mcp corpus (Postgres + pgvector)"
    assert env["data"]["target_family"] == "human kinome"
    assert env["data"]["target_count"] == 512
    assert env["data"]["chembl_release"] == "ChEMBL_34"


def test_corpus_manifest_markdown():
    md = CorpusManifest(
        target_family="human kinome",
        built_at="2026-05-30T00:00:00+00:00",
        chembl_release="ChEMBL_34",
        uniprot_snapshot="2026_02",
        protein_embedding_model="ESM-2 650M",
        chem_embedding_model=None,
        target_count=512,
        compound_count=41234,
        activity_count=87654,
    ).to_markdown()
    assert "human kinome" in md
    assert "ChEMBL_34" in md
    assert "512" in md


# ---------------------------------------------------------------------------
# Integration: real Postgres + pgvector. Runs only when GENESIS_CORPUS_DSN is set
# (the CI job provides one; local runs without Docker skip these).
# ---------------------------------------------------------------------------

_HAS_DB = bool(os.environ.get("GENESIS_CORPUS_DSN"))
_skip_no_db = pytest.mark.skipif(not _HAS_DB, reason="GENESIS_CORPUS_DSN not set")


@_skip_no_db
@pytest.mark.asyncio
async def test_integration_pool_applies_schema_on_empty_store():
    # create_corpus_pool connects + applies the schema idempotently; a fresh store has no
    # manifest row yet, so corpus_describe reports the corpus as configured-but-unbuilt.
    pool = await create_corpus_pool()
    assert pool is not None, "expected a live pool when GENESIS_CORPUS_DSN is set"
    try:
        assert await fetch_manifest(pool) is None  # empty store
        server.mcp.state = SimpleNamespace(corpus_pool=pool)
        out = await corpus_describe(CorpusDescribeInput(response_format="json"))
        env = json.loads(out)
        assert env["error"]["status"] == "NotFound"
    finally:
        await pool.close()
