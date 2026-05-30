"""CLI entry point for corpus ingestion.

    uv run python -m genesis_bio_mcp.ingest targets [LIMIT]      # kinome + ESM-2 embeddings
    uv run python -m genesis_bio_mcp.ingest activities [LIMIT]   # ChEMBL compounds + activities

Requires GENESIS_CORPUS_DSN. ``targets`` also needs the optional 'ingest' extras
(uv sync --group ingest) for ESM-2; ``activities`` needs only the core deps (RDKit). Run
``targets`` first. ``LIMIT`` caps the number of targets processed — useful for a partial build.
"""

from __future__ import annotations

import asyncio
import logging
import sys
from datetime import date

import httpx

from genesis_bio_mcp.corpus import create_corpus_pool
from genesis_bio_mcp.ingest.chembl import (
    build_compound_record,
    fetch_activities_for_target,
    resolve_chembl_target,
)
from genesis_bio_mcp.ingest.kinome import fetch_kinome_targets
from genesis_bio_mcp.ingest.load import load_activities, load_compounds, load_targets

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("genesis_bio_mcp.ingest")


async def ingest_targets(limit: int | None = None) -> int:
    pool = await create_corpus_pool()
    if pool is None:
        logger.error(
            "Corpus DB not available — set GENESIS_CORPUS_DSN and start Postgres "
            "(docker compose up -d)."
        )
        return 1
    try:
        async with httpx.AsyncClient(follow_redirects=True, timeout=60.0) as client:
            records = await fetch_kinome_targets(client, max_targets=limit)
        logger.info("Fetched %d targets; computing ESM-2 embeddings on CPU…", len(records))
        # Imported here so the heavy ML dep is only required for the embedding step.
        from genesis_bio_mcp.ingest.embed import embed_sequences

        embeddings = embed_sequences([r.sequence for r in records])
        written = await load_targets(
            pool, records, embeddings, source_version=f"UniProt {date.today().isoformat()}"
        )
        logger.info("Done — loaded %d targets into the corpus.", written)
        return 0
    finally:
        await pool.close()


async def ingest_activities(limit: int | None = None) -> int:
    pool = await create_corpus_pool()
    if pool is None:
        logger.error(
            "Corpus DB not available — set GENESIS_CORPUS_DSN and start Postgres "
            "(docker compose up -d)."
        )
        return 1
    try:
        async with pool.acquire() as conn:
            rows = await conn.fetch(
                "SELECT uniprot_accession FROM targets ORDER BY uniprot_accession"
            )
        accessions = [r["uniprot_accession"] for r in rows]
        if not accessions:
            logger.error("No targets in the corpus — run `ingest targets` first.")
            return 1
        if limit is not None:
            accessions = accessions[:limit]

        version = f"ChEMBL {date.today().isoformat()}"
        total_compounds = total_activities = 0
        async with httpx.AsyncClient(follow_redirects=True, timeout=60.0) as client:
            for acc in accessions:
                target_id = await resolve_chembl_target(client, acc)
                if not target_id:
                    logger.info("No ChEMBL target for %s — skipping.", acc)
                    continue
                activities, smiles_by_mol = await fetch_activities_for_target(
                    client, target_id, acc
                )
                if not activities:
                    continue
                # Compounds must be loaded before the activities that FK to them.
                compounds = [build_compound_record(mid, smi) for mid, smi in smiles_by_mol.items()]
                await load_compounds(pool, compounds, source_version=version)
                await load_activities(pool, activities, source_version=version)
                total_compounds += len(compounds)
                total_activities += len(activities)
                logger.info(
                    "%s (%s): %d compounds, %d activities",
                    acc,
                    target_id,
                    len(compounds),
                    len(activities),
                )
        logger.info(
            "Done — loaded %d compounds and %d activities across %d targets.",
            total_compounds,
            total_activities,
            len(accessions),
        )
        return 0
    finally:
        await pool.close()


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    command = argv[0] if argv else "targets"
    limit = int(argv[1]) if len(argv) > 1 else None
    if command == "targets":
        return asyncio.run(ingest_targets(limit))
    if command == "activities":
        return asyncio.run(ingest_activities(limit))
    logger.error(
        "Unknown command %r. Usage: python -m genesis_bio_mcp.ingest {targets|activities} [LIMIT]",
        command,
    )
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
