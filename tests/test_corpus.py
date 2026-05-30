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
from genesis_bio_mcp.server import (
    CorpusDescribeInput,
    CorpusSearchTargetsInput,
    corpus_describe,
    corpus_search_targets_by_sequence,
)


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
# corpus_search_targets_by_sequence — unit tests (DB mocked)
# ---------------------------------------------------------------------------


@pytest.mark.asyncio
async def test_search_targets_unconfigured_returns_upstream_unavailable():
    server.mcp.state = SimpleNamespace(corpus_pool=None)
    out = await corpus_search_targets_by_sequence(
        CorpusSearchTargetsInput(query="BRAF", response_format="json")
    )
    env = json.loads(out)
    assert env["error"]["status"] == "UpstreamUnavailable"


@pytest.mark.asyncio
async def test_search_targets_not_in_corpus_returns_not_found(monkeypatch):
    server.mcp.state = SimpleNamespace(corpus_pool=object())

    async def _none(_pool, _q, _k):
        return None

    monkeypatch.setattr(server, "search_similar_targets", _none)
    out = await corpus_search_targets_by_sequence(
        CorpusSearchTargetsInput(query="NOTAKINASE", response_format="json")
    )
    env = json.loads(out)
    assert env["error"]["status"] == "NotFound"


@pytest.mark.asyncio
async def test_search_targets_returns_ranked_neighbors(monkeypatch):
    server.mcp.state = SimpleNamespace(corpus_pool=object())

    async def _hits(_pool, _q, _k):
        return (
            "BRAF",
            [
                {
                    "gene_symbol": "ARAF",
                    "target_chembl_id": "CHEMBL123",
                    "uniprot_accession": "P10398",
                    "pref_name": "Serine/threonine-protein kinase A-Raf",
                    "kinase_group": "TKL",
                    "cosine_similarity": 0.99,
                    "activity_count": 42,
                }
            ],
        )

    monkeypatch.setattr(server, "search_similar_targets", _hits)
    out = await corpus_search_targets_by_sequence(
        CorpusSearchTargetsInput(query="braf", response_format="json")
    )
    env = json.loads(out)
    assert env["provenance"]["source"] == "genesis-bio-mcp corpus (Postgres + pgvector)"
    assert env["provenance"]["query"] == "BRAF"  # resolved query_gene
    assert env["data"]["hits"][0]["gene_symbol"] == "ARAF"
    assert env["data"]["hits"][0]["cosine_similarity"] == 0.99


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


def _vec(*head: float) -> list[float]:
    """A 640-d ESM-2-150M-shaped vector with the given leading components, zero-padded."""
    v = list(head) + [0.0] * (640 - len(head))
    return v


@_skip_no_db
@pytest.mark.asyncio
async def test_integration_target_similarity_knn_ranks_by_cosine():
    # Insert synthetic targets with known 640-d vectors and verify pgvector cosine kNN
    # returns the nearest neighbor first via corpus_search_targets_by_sequence.
    pool = await create_corpus_pool()
    assert pool is not None
    try:
        async with pool.acquire() as conn:
            await conn.execute("DELETE FROM targets")
            rows = [
                ("CHEMBL_BRAF", "BRAF", "P15056", _vec(1.0, 0.0)),  # query
                ("CHEMBL_ARAF", "ARAF", "P10398", _vec(0.9, 0.1)),  # near (cos ~0.993)
                ("CHEMBL_FAR", "FARGENE", "Q00000", _vec(0.0, 1.0)),  # orthogonal (cos 0)
            ]
            for tid, gene, acc, vec in rows:
                await conn.execute(
                    "INSERT INTO targets (target_chembl_id, gene_symbol, uniprot_accession, "
                    "sequence_embedding) VALUES ($1, $2, $3, $4)",
                    tid,
                    gene,
                    acc,
                    vec,
                )

        server.mcp.state = SimpleNamespace(corpus_pool=pool)
        out = await corpus_search_targets_by_sequence(
            CorpusSearchTargetsInput(query="BRAF", response_format="json")
        )
        env = json.loads(out)
        hits = env["data"]["hits"]
        assert [h["gene_symbol"] for h in hits] == ["ARAF", "FARGENE"]  # nearest first
        assert hits[0]["cosine_similarity"] > hits[1]["cosine_similarity"]
        assert hits[0]["cosine_similarity"] > 0.9
    finally:
        async with pool.acquire() as conn:
            await conn.execute("DELETE FROM targets")
        await pool.close()
