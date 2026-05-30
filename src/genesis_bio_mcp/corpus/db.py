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
       (SELECT count(*) FROM activities a WHERE a.uniprot_accession = t.uniprot_accession)
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


_COMPOUND_SQL = """
SELECT c.molecule_chembl_id,
       c.canonical_smiles,
       c.inchikey,
       c.mol_weight,
       1 - (c.morgan_fp <%> $1::bit(2048)) AS tanimoto,
       (SELECT count(*) FROM activities a WHERE a.molecule_chembl_id = c.molecule_chembl_id)
           AS activity_count
FROM compounds c
WHERE c.morgan_fp IS NOT NULL
ORDER BY c.morgan_fp <%> $1::bit(2048)
LIMIT $2
"""


async def search_similar_compounds(pool: asyncpg.Pool, fp_bits: str, top_k: int) -> list[dict]:
    """Rank corpus compounds by Tanimoto to a query Morgan fingerprint (pgvector Jaccard).

    ``fp_bits`` is a 2048-char '0'/'1' bitstring (see ``cheminformatics.morgan_fp_bits``).
    Tanimoto = 1 − Jaccard distance (``<%>``). Returns up to ``top_k`` compound dicts ordered
    by descending Tanimoto; empty if the corpus has no fingerprinted compounds.
    """
    async with pool.acquire() as conn:
        # asyncpg's bit codec needs an asyncpg.BitString, not a plain '0'/'1' str.
        rows = await conn.fetch(_COMPOUND_SQL, asyncpg.BitString(fp_bits), top_k)
    return [dict(r) for r in rows]


async def resolve_target_accession(pool: asyncpg.Pool, query: str) -> str | None:
    """Resolve a gene symbol / UniProt accession to a corpus target's accession, or None."""
    async with pool.acquire() as conn:
        return await conn.fetchval(
            "SELECT uniprot_accession FROM targets "
            "WHERE upper(gene_symbol) = upper($1) OR uniprot_accession = $1 LIMIT 1",
            query,
        )


# Hybrid retrieval: SQL filters (target / assay type / potency / MW) + optional Morgan-Tanimoto
# similarity, in ONE query. `best` picks each compound's strongest matching activity
# (DISTINCT ON … ORDER BY pchembl DESC). When a query fingerprint is supplied, results are
# ranked by pgvector Jaccard distance (most similar first); otherwise by potency.
_HYBRID_SQL = """
WITH best AS (
    SELECT DISTINCT ON (a.molecule_chembl_id)
           a.molecule_chembl_id, a.pchembl_value, a.standard_type,
           a.standard_units, a.standard_value, a.assay_confidence_score
    FROM activities a
    WHERE ($1::text IS NULL OR a.uniprot_accession = $1)
      AND ($2::text IS NULL OR a.standard_type = $2)
      AND ($3::float8 IS NULL OR a.pchembl_value >= $3)
    ORDER BY a.molecule_chembl_id, a.pchembl_value DESC NULLS LAST
)
SELECT c.molecule_chembl_id, c.canonical_smiles, c.inchikey, c.mol_weight,
       b.pchembl_value, b.standard_type, b.standard_units, b.standard_value,
       b.assay_confidence_score,
       CASE WHEN $5::bit(2048) IS NOT NULL
            THEN 1 - (c.morgan_fp <%> $5::bit(2048)) END AS tanimoto
FROM best b
JOIN compounds c ON c.molecule_chembl_id = b.molecule_chembl_id
WHERE ($4::float8 IS NULL OR c.mol_weight <= $4)
  AND ($5::bit(2048) IS NULL OR c.morgan_fp IS NOT NULL)
ORDER BY
    CASE WHEN $5::bit(2048) IS NOT NULL THEN (c.morgan_fp <%> $5::bit(2048)) END ASC NULLS LAST,
    b.pchembl_value DESC NULLS LAST
LIMIT $6 OFFSET $7
"""


async def hybrid_search_compounds(
    pool: asyncpg.Pool,
    *,
    uniprot_accession: str | None,
    standard_type: str | None,
    min_pchembl: float | None,
    max_mol_weight: float | None,
    query_fp_bits: str | None,
    limit: int,
    offset: int,
) -> list[dict]:
    """Hybrid compound search: relational filters + optional Tanimoto similarity, one query.

    Filters compounds by target / assay type / minimum pChEMBL / max MW (SQL), and — when a
    query fingerprint is given — ranks by Morgan-Tanimoto (pgvector Jaccard); otherwise by
    potency. Each compound appears once with its strongest matching activity.
    """
    fp = asyncpg.BitString(query_fp_bits) if query_fp_bits else None
    async with pool.acquire() as conn:
        rows = await conn.fetch(
            _HYBRID_SQL,
            uniprot_accession,
            standard_type,
            min_pchembl,
            max_mol_weight,
            fp,
            limit,
            offset,
        )
    return [dict(r) for r in rows]
