"""Embedding-backed corpus store (v0.6.0 hybrid retrieval).

A persistent PostgreSQL + pgvector tier over a curated bioactivity corpus (the human
kinome), queried at request time by the read-only ``corpus_*`` MCP tools. The corpus is
OPTIONAL: when ``GENESIS_CORPUS_DSN`` is unset the server runs normally and the corpus_*
tools report the store as unconfigured.

Heavy ML (ESM-2 protein embeddings) lives only in the offline ingestion tier — never on
the request path. See docs/ROADMAP.md (v0.6.0).
"""

from __future__ import annotations

from genesis_bio_mcp.corpus.db import (
    create_corpus_pool,
    fetch_manifest,
    hybrid_search_compounds,
    resolve_target_accession,
    search_similar_compounds,
    search_similar_targets,
)

__all__ = [
    "create_corpus_pool",
    "fetch_manifest",
    "hybrid_search_compounds",
    "resolve_target_accession",
    "search_similar_compounds",
    "search_similar_targets",
]
