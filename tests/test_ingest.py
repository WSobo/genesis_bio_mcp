"""Unit tests for the offline corpus ingestion pipeline (v0.6.0 M1).

The ESM-2 *model run* needs the optional `ingest` extras (torch) and isn't exercised here;
the pooling math, the UniProt fetch/parse, and the DB load are all tested without torch
(load runs against real Postgres in CI, skipped locally when GENESIS_CORPUS_DSN is unset).
"""

from __future__ import annotations

import json
import os
from types import SimpleNamespace

import httpx
import numpy as np
import pytest
import respx

from genesis_bio_mcp import server
from genesis_bio_mcp.corpus import create_corpus_pool
from genesis_bio_mcp.ingest.embed import EMBED_DIM, _mean_pool, embed_sequences
from genesis_bio_mcp.ingest.kinome import (
    TargetRecord,
    _next_link,
    _parse_entry,
    fetch_kinome_targets,
)
from genesis_bio_mcp.ingest.load import load_targets
from genesis_bio_mcp.server import CorpusDescribeInput, corpus_describe


def test_mean_pool_excludes_bos_and_eos():
    # 4 token rows: [BOS, res1, res2, EOS]; mean-pool should average res1 + res2 only.
    reps = np.array([[9.0, 9.0], [1.0, 2.0], [3.0, 4.0], [9.0, 9.0]], dtype="float32")
    pooled = _mean_pool(reps, seq_len=2)
    assert pooled == [2.0, 3.0]  # mean of [1,2] and [3,4]


def test_mean_pool_degenerate_falls_back_to_all_rows():
    reps = np.array([[1.0, 1.0]], dtype="float32")
    assert _mean_pool(reps, seq_len=0) == [1.0, 1.0]


def test_embed_sequences_requires_ingest_extras():
    # torch / fair-esm are not in the default env → actionable ImportError, not a bare crash.
    with pytest.raises(ImportError, match="uv sync --group ingest"):
        embed_sequences(["MKKFFDSR"])


def test_parse_entry_extracts_fields():
    entry = {
        "primaryAccession": "P15056",
        "sequence": {"value": "MAALSGGGGGG"},
        "genes": [{"geneName": {"value": "BRAF"}}],
        "proteinDescription": {"recommendedName": {"fullName": {"value": "Serine/threonine BRAF"}}},
        "organism": {"scientificName": "Homo sapiens"},
    }
    rec = _parse_entry(entry)
    assert rec is not None
    assert rec.uniprot_accession == "P15056"
    assert rec.gene_symbol == "BRAF"
    assert rec.sequence == "MAALSGGGGGG"


def test_parse_entry_skips_entries_without_sequence():
    assert _parse_entry({"primaryAccession": "X", "genes": []}) is None


def test_next_link_parses_cursor():
    header = '<https://rest.uniprot.org/uniprotkb/search?cursor=abc&size=500>; rel="next"'
    assert _next_link(header) == "https://rest.uniprot.org/uniprotkb/search?cursor=abc&size=500"
    assert _next_link(None) is None
    assert _next_link('<x>; rel="last"') is None


@pytest.mark.asyncio
@respx.mock
async def test_fetch_kinome_targets_single_page():
    payload = {
        "results": [
            {
                "primaryAccession": "P15056",
                "sequence": {"value": "MAAL"},
                "genes": [{"geneName": {"value": "BRAF"}}],
                "proteinDescription": {"recommendedName": {"fullName": {"value": "BRAF"}}},
                "organism": {"scientificName": "Homo sapiens"},
            },
            {"primaryAccession": "NOSEQ", "genes": []},  # dropped — no sequence
        ]
    }
    respx.get(url__startswith="https://rest.uniprot.org/uniprotkb/search").mock(
        return_value=httpx.Response(200, json=payload)  # no Link header → single page
    )
    async with httpx.AsyncClient() as client:
        records = await fetch_kinome_targets(client)
    assert [r.gene_symbol for r in records] == ["BRAF"]


# ---------------------------------------------------------------------------
# Integration: load targets into real Postgres + pgvector (CI only).
# ---------------------------------------------------------------------------

_skip_no_db = pytest.mark.skipif(
    not os.environ.get("GENESIS_CORPUS_DSN"), reason="GENESIS_CORPUS_DSN not set"
)


@_skip_no_db
@pytest.mark.asyncio
async def test_integration_load_targets_and_manifest():
    pool = await create_corpus_pool()
    assert pool is not None
    try:
        async with pool.acquire() as conn:
            await conn.execute("DELETE FROM targets")
        records = [
            TargetRecord("P15056", "BRAF", "Serine/threonine BRAF", "Homo sapiens", "MAAL"),
            TargetRecord("P10398", "ARAF", "A-Raf", "Homo sapiens", "MGNG"),
        ]
        embeddings = [[0.1] * EMBED_DIM, [0.2] * EMBED_DIM]
        written = await load_targets(pool, records, embeddings, source_version="UniProt test")
        assert written == 2

        # The manifest now reports the loaded targets via corpus_describe.
        server.mcp.state = SimpleNamespace(corpus_pool=pool)
        out = await corpus_describe(CorpusDescribeInput(response_format="json"))
        env = json.loads(out)
        assert env["data"]["target_count"] == 2
        assert env["data"]["target_family"] == "human kinome"
        assert "ESM-2-150M" in env["data"]["protein_embedding_model"]
    finally:
        async with pool.acquire() as conn:
            await conn.execute("DELETE FROM targets")
        await pool.close()
