"""Async Postgres + pgvector access for the corpus store.

The corpus is OPTIONAL. ``create_corpus_pool`` returns ``None`` when no DSN is configured
or the database is unreachable — it NEVER raises — so the rest of the server always starts.
Schema application is idempotent (see ``schema.sql``).
"""

from __future__ import annotations

import logging
from importlib import resources

import asyncpg

from genesis_bio_mcp.config.settings import settings

logger = logging.getLogger(__name__)


def _schema_sql() -> str:
    """Read the bundled corpus DDL."""
    return (
        resources.files("genesis_bio_mcp.corpus").joinpath("schema.sql").read_text(encoding="utf-8")
    )


async def create_corpus_pool() -> asyncpg.Pool | None:
    """Create the asyncpg pool and ensure the schema, or return ``None``.

    Returns ``None`` (never raises) when ``GENESIS_CORPUS_DSN`` is unset, or when a DSN is
    configured but the database cannot be reached / the schema cannot be applied — logging a
    warning in the latter case. This keeps the corpus tier strictly optional.
    """
    dsn = settings.corpus_dsn
    if not dsn:
        logger.info(
            "GENESIS_CORPUS_DSN not set — corpus_* tools will report the store as unconfigured."
        )
        return None
    try:
        pool = await asyncpg.create_pool(
            dsn,
            min_size=settings.corpus_pool_min_size,
            max_size=settings.corpus_pool_max_size,
        )
        async with pool.acquire() as conn:
            await conn.execute(_schema_sql())
        logger.info("Corpus store connected; schema ensured.")
        return pool
    except Exception as exc:
        logger.warning(
            "Corpus store unavailable (GENESIS_CORPUS_DSN is set but connect/schema failed): %s",
            exc,
        )
        return None


async def fetch_manifest(pool: asyncpg.Pool) -> dict | None:
    """Return the single corpus-manifest row as a dict, or ``None`` if the corpus is unbuilt."""
    row = await pool.fetchrow("SELECT * FROM corpus_manifest WHERE id = 1")
    return dict(row) if row else None
