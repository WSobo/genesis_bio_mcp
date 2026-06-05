"""Regression tests: clients must NOT cache transient failures (cache poisoning).

Before these fixes, gnomAD/MaveDB cached ``None``/``[]`` inside their ``except``
blocks and MyVariant/UniProt cached ``None`` unconditionally — so a single network
blip permanently masked that gene/URN for the rest of the server session. Each test
forces the HTTP layer to raise, then asserts the session cache was left empty so the
next call can re-fetch.
"""

from __future__ import annotations

from unittest.mock import AsyncMock

import httpx

from genesis_bio_mcp.clients.gnomad import GnomADClient
from genesis_bio_mcp.clients.mavedb import MaveDBClient
from genesis_bio_mcp.clients.myvariant import MyVariantClient
from genesis_bio_mcp.clients.uniprot import UniProtClient


def _failing_client() -> httpx.AsyncClient:
    """A real AsyncClient whose get/post always raise a transient network error."""
    client = httpx.AsyncClient()
    client.get = AsyncMock(side_effect=httpx.ConnectError("boom"))  # type: ignore[method-assign]
    client.post = AsyncMock(side_effect=httpx.ConnectError("boom"))  # type: ignore[method-assign]
    return client


async def test_myvariant_get_annotation_does_not_cache_failure() -> None:
    http = _failing_client()
    try:
        mv = MyVariantClient(http)
        result = await mv.get_annotation("chr17:g.7676154G>A")
        assert result is None
        assert mv._cache == {}  # not poisoned — a retry is still possible
    finally:
        await http.aclose()


async def test_uniprot_get_protein_does_not_cache_failure() -> None:
    http = _failing_client()
    try:
        uni = UniProtClient(http)
        result = await uni.get_protein("EGFR")
        assert result is None
        assert uni._cache == {}
    finally:
        await http.aclose()


async def test_gnomad_load_variants_does_not_cache_failure() -> None:
    http = _failing_client()
    try:
        gn = GnomADClient(http)
        result = await gn._load_variants("TP53")
        assert result is None
        assert gn._variants_cache == {}
    finally:
        await http.aclose()


async def test_mavedb_load_scores_does_not_cache_failure() -> None:
    http = _failing_client()
    try:
        mave = MaveDBClient(http)
        result = await mave._load_scores("urn:mavedb:00000001-a-1")
        assert result == []
        assert mave._scores_cache == {}
    finally:
        await http.aclose()
