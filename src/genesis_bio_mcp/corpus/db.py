"""Async Postgres + pgvector access for the corpus store.

The corpus is OPTIONAL. ``create_corpus_pool`` returns ``None`` when no DSN is configured
or the database is unreachable — it NEVER raises — so the rest of the server always starts.
Schema application is idempotent (see ``schema.sql``).
"""

from __future__ import annotations

import logging
from importlib import resources

import asyncpg
from pgvector.asyncpg import register_vector

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
        # Ensure the pgvector extension + schema on a one-off connection FIRST, so the
        # 'vector' type exists before the pool registers its codec on every connection.
        conn = await asyncpg.connect(dsn)
        try:
            await conn.execute(_schema_sql())
        finally:
            await conn.close()

        async def _init(connection: asyncpg.Connection) -> None:
            # Register the pgvector codec so vector columns round-trip as Python lists /
            # numpy arrays and can be passed straight back into similarity queries.
            await register_vector(connection)

        pool = await asyncpg.create_pool(
            dsn,
            min_size=settings.corpus_pool_min_size,
            max_size=settings.corpus_pool_max_size,
            init=_init,
        )
        logger.info("Corpus store connected; schema ensured; pgvector codec registered.")
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


_NEIGHBOR_SQL = """
SELECT t.target_chembl_id,
       t.gene_symbol,
       t.uniprot_accession,
       t.pref_name,
       t.kinase_group,
       1 - (t.sequence_embedding <=> $1) AS cosine_similarity,
       (SELECT count(*) FROM activities a WHERE a.target_chembl_id = t.target_chembl_id)
           AS activity_count
FROM targets t
WHERE t.gene_symbol <> $2 AND t.sequence_embedding IS NOT NULL
ORDER BY t.sequence_embedding <=> $1
LIMIT $3
"""


async def search_similar_targets(
    pool: asyncpg.Pool, query: str, top_k: int
) -> tuple[str, list[dict]] | None:
    """Find targets whose ESM-2 embedding is closest to the query target's STORED vector.

    The query (gene symbol or UniProt accession) is resolved to a vector already in the
    corpus — no embedding happens at request time. Ranks the rest by pgvector cosine
    distance (``<=>``). Returns ``(resolved_gene_symbol, hits)``, or ``None`` if the query
    target is not in the corpus / has no embedding.
    """
    async with pool.acquire() as conn:
        row = await conn.fetchrow(
            "SELECT gene_symbol, sequence_embedding FROM targets "
            "WHERE (upper(gene_symbol) = upper($1) OR uniprot_accession = $1) "
            "AND sequence_embedding IS NOT NULL LIMIT 1",
            query,
        )
        if row is None:
            return None
        hits = await conn.fetch(_NEIGHBOR_SQL, row["sequence_embedding"], row["gene_symbol"], top_k)
    return row["gene_symbol"], [dict(h) for h in hits]
