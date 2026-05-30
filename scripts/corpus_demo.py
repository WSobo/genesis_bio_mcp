"""End-to-end demo of the embedding-backed corpus (`corpus_*`) tools.

Requires a built corpus and GENESIS_CORPUS_DSN set (see docs/ROADMAP.md v0.6.0 /
docs/corpus-eval.md). Read-only; prints each tool's Markdown to stdout.

    docker compose up -d
    export GENESIS_CORPUS_DSN="postgresql://genesis:genesis@localhost:5432/genesis_corpus"
    uv run python -m genesis_bio_mcp.ingest targets 50 && uv run python -m genesis_bio_mcp.ingest activities 50
    uv run python scripts/corpus_demo.py BRAF "CC(=O)Oc1ccccc1C(=O)O"

Args: [GENE] [SMILES] — default BRAF and aspirin.
"""

from __future__ import annotations

import asyncio
import sys
from types import SimpleNamespace

from genesis_bio_mcp import server
from genesis_bio_mcp.corpus import create_corpus_pool
from genesis_bio_mcp.server import (
    CorpusDescribeInput,
    CorpusFindSimilarCompoundsInput,
    CorpusSearchCompoundsInput,
    CorpusSearchTargetsInput,
    corpus_describe,
    corpus_find_similar_compounds,
    corpus_search_compounds,
    corpus_search_targets_by_sequence,
)


def _banner(title: str) -> None:
    print(f"\n{'=' * 78}\n{title}\n{'=' * 78}")


async def main(gene: str, smiles: str) -> int:
    pool = await create_corpus_pool()
    if pool is None:
        print("Corpus not available — set GENESIS_CORPUS_DSN and `docker compose up -d`.")
        return 1
    server.mcp.state = SimpleNamespace(corpus_pool=pool)
    try:
        _banner("corpus_describe — coverage + provenance")
        print(await corpus_describe(CorpusDescribeInput()))

        _banner(f"corpus_search_targets_by_sequence — targets similar to {gene}")
        print(await corpus_search_targets_by_sequence(CorpusSearchTargetsInput(query=gene)))

        _banner(f"corpus_find_similar_compounds — compounds like {smiles}")
        print(await corpus_find_similar_compounds(CorpusFindSimilarCompoundsInput(smiles=smiles)))

        _banner(f"corpus_search_compounds — potent compounds vs {gene} (pChEMBL ≥ 7)")
        print(
            await corpus_search_compounds(
                CorpusSearchCompoundsInput(target=gene, min_pchembl=7.0, limit=10)
            )
        )

        _banner(f"corpus_search_compounds — HYBRID: vs {gene}, ranked by similarity to query")
        print(
            await corpus_search_compounds(
                CorpusSearchCompoundsInput(target=gene, similar_to_smiles=smiles, limit=10)
            )
        )
        return 0
    finally:
        await pool.close()


if __name__ == "__main__":
    gene_arg = sys.argv[1] if len(sys.argv) > 1 else "BRAF"
    smiles_arg = sys.argv[2] if len(sys.argv) > 2 else "CC(=O)Oc1ccccc1C(=O)O"
    raise SystemExit(asyncio.run(main(gene_arg, smiles_arg)))
