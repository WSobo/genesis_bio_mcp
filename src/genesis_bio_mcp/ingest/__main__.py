"""CLI entry point for corpus ingestion.

    uv run python -m genesis_bio_mcp.ingest targets [LIMIT]

Requires GENESIS_CORPUS_DSN and the optional 'ingest' extras (uv sync --group ingest).
``LIMIT`` (optional) caps the number of targets — useful for a quick partial build.
"""

from __future__ import annotations

import asyncio
import logging
import sys
from datetime import date

import httpx

from genesis_bio_mcp.corpus import create_corpus_pool
from genesis_bio_mcp.ingest.kinome import fetch_kinome_targets
from genesis_bio_mcp.ingest.load import load_targets

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


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    command = argv[0] if argv else "targets"
    limit = int(argv[1]) if len(argv) > 1 else None
    if command == "targets":
        return asyncio.run(ingest_targets(limit))
    logger.error(
        "Unknown command %r. Usage: python -m genesis_bio_mcp.ingest targets [LIMIT]", command
    )
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
