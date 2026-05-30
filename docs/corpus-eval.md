# Corpus evaluation suite (v0.6.0)

Ten evaluation question/answer pairs for the embedding-backed corpus (`corpus_*`) tools,
following MCP eval best practice — each is **independent** (no ordering dependency),
**read-only**, **multi-call or filter-bearing**, and **verifiable** (a checkable expectation,
not a vibe).

These run against a **built corpus** — bring one up first (a small partial build is enough):

```bash
docker compose up -d
export GENESIS_CORPUS_DSN="postgresql://genesis:genesis@localhost:5432/genesis_corpus"
uv sync --group ingest
uv run python -m genesis_bio_mcp.ingest targets 50       # ~50 kinases + ESM-2 embeddings
uv run python -m genesis_bio_mcp.ingest activities 50    # their ChEMBL compounds + activities
```

Expectations below are phrased against the corpus contents, so exact gene/compound hits depend
on which slice you ingested; the **verifiable property** (filter obeyed, ordering, typed error)
holds regardless.

| # | Question | Tool call(s) | Verifiable expectation |
|---|---|---|---|
| 1 | What does the corpus cover and how was it built? | `corpus_describe` | `target_family` = "human kinome"; non-zero target/compound/activity counts; `protein_embedding_model` mentions ESM-2-150M; a `built_at` timestamp + ChEMBL/UniProt provenance. |
| 2 | Which targets are most similar to BRAF by protein sequence? | `corpus_search_targets_by_sequence(query="BRAF")` | Returns neighbors ranked by descending cosine similarity; BRAF's RAF-family / TKL-group relatives (e.g. ARAF, RAF1) rank near the top; every hit has `cosine_similarity` ≤ 1. |
| 3 | What chemical matter in the corpus is structurally like aspirin? | `corpus_find_similar_compounds(smiles="CC(=O)Oc1ccccc1C(=O)O")` | Hits ranked by descending Tanimoto, each in [0,1]; if aspirin itself is indexed it appears at Tanimoto ≈ 1.0. |
| 4 | What are the most potent compounds against a given kinase? | `corpus_search_compounds(target="<gene>")` | Default ranking is by descending pChEMBL; the top hit has the highest pChEMBL of the returned set. |
| 5 | Find sub-µM IC50 compounds against a kinase. | `corpus_search_compounds(target="<gene>", standard_type="IC50", min_pchembl=6)` | **Every** returned hit has `standard_type` = IC50 and `pchembl_value` ≥ 6. |
| 6 | Find lead-like (MW ≤ 500), potent (pChEMBL ≥ 7) compounds for a target. | `corpus_search_compounds(target="<gene>", max_mol_weight=500, min_pchembl=7)` | Every hit has `mol_weight` ≤ 500 and `pchembl_value` ≥ 7. |
| 7 | Compounds structurally similar to a query that *also* have measured activity against a target (the hybrid question). | `corpus_search_compounds(target="<gene>", similar_to_smiles="<SMILES>")` | Results ranked by descending Tanimoto; every hit also has a `pchembl_value` against the target (similarity AND activity, together). |
| 8 | Is a gene covered, and if so what inhibitors exist? (multi-call) | `corpus_describe` → `corpus_search_compounds(target="<gene>")` | If the gene is in the corpus, the compound search returns bioactive compounds; the manifest counts are consistent with non-empty results. |
| 9 | For a query target, find its nearest corpus neighbor and that neighbor's chemical matter. (multi-call) | `corpus_search_targets_by_sequence(query="<gene>")` → `corpus_search_compounds(target="<neighbor>")` | The second call's `target` is a gene returned by the first; it yields that neighbor's bioactive compounds. |
| 10 | Search compounds for a target that is **not** in the corpus (robustness). | `corpus_search_compounds(target="ALB")` (albumin — not a kinase) | JSON error envelope with `error.status` = `"NotFound"` and a message pointing to `corpus_describe`; no crash. |

**Notes for graders**
- Use `response_format="json"` to assert on structured fields (`provenance`, `data`, `error.status`).
- Tools degrade deterministically: no `GENESIS_CORPUS_DSN` → `error.status` = `UpstreamUnavailable`;
  corpus configured-but-unbuilt → `NotFound`. These are themselves verifiable behaviors.
- Similarity (Tanimoto / cosine) is a **retrieval** signal, not a validated activity or functional
  prediction — graders should not score it as a binding claim.

A runnable end-to-end demonstration of tools 1–7 is in
[`scripts/corpus_demo.py`](../scripts/corpus_demo.py).
