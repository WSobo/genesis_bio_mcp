"""Load target rows + embeddings into Postgres and refresh the corpus manifest.

Pure asyncpg writes — no ML here. The pgvector codec is registered on the pool
(see ``corpus.db.create_corpus_pool``), so embeddings can be passed as plain lists.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from datetime import UTC, datetime

import asyncpg

from genesis_bio_mcp.ingest.embed import ESM2_MODEL
from genesis_bio_mcp.ingest.kinome import TargetRecord

logger = logging.getLogger(__name__)


@dataclass
class CompoundRecord:
    """A compound to index: ChEMBL-keyed, standardized, with its Morgan fingerprint bits."""

    molecule_chembl_id: str
    canonical_smiles: str | None
    inchi: str | None
    inchikey: str | None
    mol_weight: float | None
    morgan_fp_bits: str | None  # 2048-char '0'/'1' bitstring, or None if unparseable


_UPSERT_TARGET = """
INSERT INTO targets (uniprot_accession, gene_symbol, pref_name, organism, sequence,
                     sequence_embedding, source, source_version, retrieved_at)
VALUES ($1, $2, $3, $4, $5, $6, $7, $8, $9)
ON CONFLICT (uniprot_accession) DO UPDATE SET
    gene_symbol = EXCLUDED.gene_symbol,
    pref_name = EXCLUDED.pref_name,
    organism = EXCLUDED.organism,
    sequence = EXCLUDED.sequence,
    sequence_embedding = EXCLUDED.sequence_embedding,
    source = EXCLUDED.source,
    source_version = EXCLUDED.source_version,
    retrieved_at = EXCLUDED.retrieved_at
"""

_UPSERT_MANIFEST = """
INSERT INTO corpus_manifest (id, built_at, target_family, protein_embedding_model,
                             target_count, compound_count, activity_count)
VALUES (1, $1, $2, $3, $4, $5, $6)
ON CONFLICT (id) DO UPDATE SET
    built_at = EXCLUDED.built_at,
    target_family = EXCLUDED.target_family,
    protein_embedding_model = EXCLUDED.protein_embedding_model,
    target_count = EXCLUDED.target_count,
    compound_count = EXCLUDED.compound_count,
    activity_count = EXCLUDED.activity_count
"""


async def _refresh_manifest(conn: asyncpg.Connection, now: datetime, target_family: str) -> None:
    counts = await conn.fetchrow(
        "SELECT (SELECT count(*) FROM targets) AS t, "
        "(SELECT count(*) FROM compounds) AS c, "
        "(SELECT count(*) FROM activities) AS a"
    )
    await conn.execute(
        _UPSERT_MANIFEST,
        now,
        target_family,
        f"ESM-2-150M ({ESM2_MODEL}), mean-pooled",
        counts["t"],
        counts["c"],
        counts["a"],
    )


async def load_targets(
    pool: asyncpg.Pool,
    records: list[TargetRecord],
    embeddings: list[list[float]],
    *,
    source_version: str,
    target_family: str = "human kinome",
) -> int:
    """Upsert targets (with embeddings) in one transaction, then refresh the manifest.

    Returns the number of target rows written. ``records`` and ``embeddings`` must align.
    """
    if len(records) != len(embeddings):
        raise ValueError(f"records ({len(records)}) and embeddings ({len(embeddings)}) misaligned")
    now = datetime.now(UTC)
    async with pool.acquire() as conn, conn.transaction():
        for rec, emb in zip(records, embeddings):
            await conn.execute(
                _UPSERT_TARGET,
                rec.uniprot_accession,
                rec.gene_symbol,
                rec.protein_name,
                rec.organism,
                rec.sequence,
                emb,
                "UniProt",
                source_version,
                now,
            )
        await _refresh_manifest(conn, now, target_family)
    logger.info("Loaded %d targets into the corpus.", len(records))
    return len(records)


_UPSERT_COMPOUND = """
INSERT INTO compounds (molecule_chembl_id, canonical_smiles, inchi, inchikey, mol_weight,
                       morgan_fp, source, source_version, retrieved_at)
VALUES ($1, $2, $3, $4, $5, $6::bit(2048), $7, $8, $9)
ON CONFLICT (molecule_chembl_id) DO UPDATE SET
    canonical_smiles = EXCLUDED.canonical_smiles,
    inchi = EXCLUDED.inchi,
    inchikey = EXCLUDED.inchikey,
    mol_weight = EXCLUDED.mol_weight,
    morgan_fp = EXCLUDED.morgan_fp,
    source = EXCLUDED.source,
    source_version = EXCLUDED.source_version,
    retrieved_at = EXCLUDED.retrieved_at
"""


async def load_compounds(
    pool: asyncpg.Pool,
    records: list[CompoundRecord],
    *,
    source_version: str,
    target_family: str = "human kinome",
) -> int:
    """Upsert compound rows (with Morgan fingerprints) in one transaction, refresh the manifest.

    The fingerprint is stored as a ``bit(2048)`` (passed as a '0'/'1' bitstring, cast in SQL),
    so pgvector Jaccard (= 1 − Tanimoto) similarity search works directly on it. Returns the
    number of compound rows written.
    """
    now = datetime.now(UTC)
    async with pool.acquire() as conn, conn.transaction():
        for rec in records:
            await conn.execute(
                _UPSERT_COMPOUND,
                rec.molecule_chembl_id,
                rec.canonical_smiles,
                rec.inchi,
                rec.inchikey,
                rec.mol_weight,
                rec.morgan_fp_bits,
                "ChEMBL",
                source_version,
                now,
            )
        await _refresh_manifest(conn, now, target_family)
    logger.info("Loaded %d compounds into the corpus.", len(records))
    return len(records)
