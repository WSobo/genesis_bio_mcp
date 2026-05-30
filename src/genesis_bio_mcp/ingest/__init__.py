"""Offline corpus ingestion (v0.6.0).

Extract → Transform → Load pipeline that builds the embedding-backed corpus into
Postgres + pgvector. Runs *off* the request path: heavy ML (ESM-2 via torch/fair-esm)
lives here only, installed via the optional ``ingest`` dependency group
(``uv sync --group ingest``) and imported lazily so the serving install never loads it.

Build the targets corpus (CPU):

    uv sync --group ingest
    export GENESIS_CORPUS_DSN="postgresql://genesis:genesis@localhost:5432/genesis_corpus"
    uv run python -m genesis_bio_mcp.ingest targets
"""
