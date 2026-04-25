# Architecture

This document covers the internal structure of genesis-bio-mcp: directory
layout, client/model patterns, design rationale, the `prioritize_target`
scoring model, and per-database API notes. For external-facing tool
documentation see [tools.md](tools.md). For running / deploying the server
see [deployment.md](deployment.md).

---

## Directory layout

```
src/genesis_bio_mcp/
├── server.py                       # FastMCP server: 27 tools + run_biology_workflow,
│                                     lifespan wiring, tool://registry resource
├── models.py                       # Pydantic V2 output models — every one implements
│                                     to_markdown() + model_dump_json()
├── workflow_agent.py               # ToolSpec registry, run_agent_loop(),
│                                     format_registry_docs() for tool://registry
├── config/
│   ├── settings.py                 # Pydantic Settings — GENESIS_* env var loader
│   ├── efo_resolver.py             # OLS4 EFO ontology client (free-text → EFO terms)
│   └── trait_synonyms.py           # Fallback GWAS trait synonym table
├── clients/                        # One file per external API
│   ├── uniprot.py                  # UniProt REST: annotation + FASTA + DISULFID features
│   ├── open_targets.py             # Open Targets v4 GraphQL: 3-step target/disease resolution
│   ├── depmap.py                   # DepMap task API + 7-day disk cache + OT proxy fallback
│   ├── gwas.py                     # GWAS Catalog HAL/REST: concurrent fetch paths + disk cache
│   ├── pubchem.py                  # PubChem REST: tenacity retries + Semaphore
│   ├── chembl.py                   # ChEMBL REST: target lookup → pChEMBL bioactivities
│   ├── alphafold.py                # AlphaFold + RCSB PDB
│   ├── string_db.py                # STRING network (score ≥700)
│   ├── biogrid.py                  # BioGRID REST (auth: BIOGRID_ACCESS_KEY)
│   ├── dgidb.py                    # DGIdb 5.0 GraphQL
│   ├── clinical_trials.py          # ClinicalTrials.gov v2 via curl_cffi (Chrome TLS impersonation)
│   ├── reactome.py                 # Reactome pathway API
│   ├── sabdab.py                   # SAbDab TSV cache + RCSB FASTA + AbNum CDR annotation
│   ├── iedb.py                     # IEDB B-cell epitope database
│   ├── iedb_tools.py               # IEDB NextGen Tools MHC-I/II binding prediction
│   ├── gnomad.py                   # gnomAD v4 GraphQL (constraint + variant lookup)
│   ├── interpro.py                 # InterPro domain annotations
│   ├── mavedb.py                   # MaveDB score-set metadata + per-variant scores
│   ├── myvariant.py                # MyVariant.info: ClinVar + AlphaMissense + REVEL + CADD aggregator
│   ├── ensembl.py                  # Ensembl REST: HGNC → ENSG, region overlap, VEP HGVS consequences
│   ├── gtex.py                     # GTEx v8: bulk-RNA median TPM per tissue (uses shared EnsemblClient)
│   ├── hpa.py                      # Human Protein Atlas: tissue-specificity, subcellular, pathology (7-day disk cache)
│   ├── openfda.py                  # OpenFDA: FAERS adverse events + boxed warnings + recalls (optional OPENFDA_API_KEY)
│   └── variant_effects.py          # Aggregator combining gnomAD + MyVariant + MaveDB + Ensembl VEP per-mutation
└── tools/
    ├── gene_resolver.py            # Multi-source alias resolution
    ├── target_prioritization.py    # asyncio.gather orchestration + scoring + confidence CI
    ├── biochem.py                  # Pure-Python MW/pI/GRAVY/ε₂₈₀ + liability-motif scanner
    └── variant_parser.py           # Parse R175H / p.Arg175His / Arg175His → canonical form
```

---

## Client pattern

Every API client follows the same contract so that extending the MCP is a
pattern-match exercise rather than a design problem.

```python
class MyApiClient:
    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client                       # shared AsyncClient — no per-client client
        self._cache: dict[str, MyModel | None] = {} # session cache keyed by input

    async def get_data(self, gene_symbol: str) -> MyModel | None:
        if gene_symbol in self._cache:
            return self._cache[gene_symbol]
        async with _SEMAPHORE:                      # rate-limit guard per API
            try:
                resp = await self._client.get(URL, timeout=20.0)
                resp.raise_for_status()
                result = _parse(resp.json())
            except Exception as exc:
                logger.warning("MyApi failed for %s: %s", gene_symbol, exc)
                return None                         # never raise to the caller
        self._cache[gene_symbol] = result
        return result
```

Non-negotiable rules: accept a shared `httpx.AsyncClient`, module-level
`asyncio.Semaphore` with a documented limit, session cache on the instance,
every HTTP call wrapped in `try/except`, explicit `timeout=`, return
`None` / `[]` on error with a `logger.warning`.

API-key-guarded clients check `os.environ.get(...)` at the top of the
public method and return `None` + warn if absent — tools degrade
gracefully rather than crashing the server.

---

## Model pattern

Every output model inherits `pydantic.BaseModel` and implements
`to_markdown() -> str`. Tools always return the Markdown form by default and
the `model_dump_json()` form when `response_format="json"` — never raw
Python dicts.

See [tools.md](tools.md) for which model each tool produces.

---

## Design decisions

| Decision | Rationale |
|----------|-----------|
| Single shared `httpx.AsyncClient` via lifespan | Connection pooling across every tool; no per-request TLS handshake overhead |
| `asyncio.gather` for parallel sub-queries | `prioritize_target` hits 6+ APIs simultaneously — latency ≈ max(single calls), not sum |
| `_safe()` / try-except on every coroutine | One API failure never crashes the pipeline; errors surface as `data_gaps` / `errors` |
| EFO ontology-backed GWAS trait matching | Free-text queries (`"fat"`, `"joint inflammation"`) resolve via OLS4 to canonical EFO terms. Covers the full GWAS vocabulary without manual curation |
| GWAS concurrent fetch paths + session gene cache | Primary (gene-ID) and SNP paths run via `asyncio.wait(timeout=15s)`. Same gene queried for multiple traits fetches associations once and filters per trait (COX2→PTGS2 pain reuses PTGS2 inflammation: 41.7s → 1.9s) |
| GWAS score saturation at 3 hits | API pagination means hit counts can shift between runs. Saturating stabilizes the scoring signal |
| OT **biologics floor** | Pharmacologically-validated biologics targets (HER2/breast, TNF/RA) would be underscored because OT's multi-datatype composite is depressed by irrelevant evidence classes. Floor at 3.25 applies only when `known_drug_score > 0.9` **and** both genetic/somatic scores are null |
| Pydantic V2 + `to_markdown()` | Tools output agent-readable Markdown; MCP clients render directly without post-processing |
| DepMap + EFO 7-day disk caches | EFO IDs are stable across quarterly releases; DepMap's ~10 MB CSV download is ~30s — cached so warm starts are instant |
| `asyncio.Semaphore` per rate-limited API | STRING / ChEMBL / BioGRID: `Semaphore(2)`. Reactome / PubChem / gnomAD: `Semaphore(3)`. Prevents 429s without serializing the pipeline |
| `response_format` on every tool | Markdown for agents, JSON for pipelines — no separate tool variants needed |
| `ToolSpec` registry with `use_when` | Every tool has an embedding-searchable sentence so `run_biology_workflow` can do semantic retrieval at scale rather than prompting Claude with all 27 descriptions |
| Shared `tools/biochem.py` | MW / pI / GRAVY / ε₂₈₀ / liability scanner used by `get_protein_sequence` and (planned) antibody developability. Avoids duplication when M4 ships |
| `curl_cffi` for ClinicalTrials.gov only | Cloudflare TLS-fingerprints stock httpx from WSL2 and some cloud IP ranges. Narrowly scoped to this one client — everything else still uses the shared httpx |

---

## Scoring model — `prioritize_target`

The composite priority score (0–10) combines seven evidence axes:

| Source | Max | Logic |
|--------|-----|-------|
| Open Targets association | 3.0† | `overall_score × 3` |
| DepMap CRISPR dependency | 2.0 | `fraction_dependent × 2` (×0.7 if OT proxy; ×1.2 if indication matches top lineage) |
| GWAS evidence | 2.0 | `min(hits, 3) / 3 × 2` — saturates at 3 replicated hits |
| Clinical / known-drug evidence | 1.5 | `known_drug_score × 1.5` |
| ChEMBL potency | 1.5 | pChEMBL ≥9 → 1.5, ≥7 → 1.0, ≥5 → 0.5, else 0.25 |
| UniProt protein quality | 1.5 | reviewed (+0.5) + variant coverage (max +1.0) |
| HPA tissue specificity (extended mode) | 1.0 | `Tissue enriched` 1.0 → `Group enriched` 0.7 → `Tissue enhanced` 0.5 → `Low tissue specificity` 0.2 → `Not detected` 0.0. Sourced from `get_protein_atlas`; only contributes when extended-mode fan-out fetches HPA |

The seven per-axis maxes sum to 12.5; the final priority score is capped at
10.0. The expression axis was added in v0.3.0 — by design it does not push
existing high-tier targets over the cap, so v0.2.x scores remain backward
compatible.

**Pan-essential cap:** pan-essential genes (DepMap `common_essential`) have
their DepMap contribution capped at 0.5 — a narrow therapeutic window is a
liability, not an asset.

†**Biologics floor:** when `known_drug_score > 0.9` **and** both
`genetic_association_score` and `somatic_mutation_score` are null, the OT
contribution is floored at 3.25. This prevents pharmacologically-validated
biologics targets (HER2/breast cancer, TNF/RA) from being underscored
because OT's multi-datatype composite is depressed by evidence classes
irrelevant to their mechanism. Targets with populated genetic or somatic
scores (BRAF/melanoma) are unaffected.

All scoring constants are documented inline in
`src/genesis_bio_mcp/tools/target_prioritization.py` with their tuning
ranges. Calibrated against the [benchmark matrix](benchmark.md) and
intentionally conservative to avoid overfitting to any single target class.

---

## API reference

| Database | API Type | Auth | Rate / Notes |
|----------|----------|------|--------------|
| UniProt | REST | — | Generous. `organism_id:9606 AND reviewed:true` for human Swiss-Prot. FASTA endpoint separate |
| Open Targets v4 | GraphQL | — | ~2 req/s. 3-step Ensembl ID + EFO ID resolution (automatic) |
| DepMap | REST (task queue) | — | Moderate. Celery polling → pre-signed CSV URL. 7-day disk cache. OT fallback if unreachable |
| GWAS Catalog | REST / HAL | — | Moderate. Concurrent gene-ID + SNP fetch paths; 15s global timeout; 24h disk cache |
| EBI OLS4 (EFO) | REST | — | Generous. 7-day cache. Ontology data queried at runtime — not redistributed |
| PubChem | REST | — | 5 req/s. HTTP 503 on rate limit; tenacity retries + `Semaphore(3)` |
| ChEMBL | REST | — | ~1 req/s. Two-step: target lookup → bioactivity. `Semaphore(2)` |
| AlphaFold | REST | — | Generous. pLDDT ≥90 = high confidence; <70 = disordered |
| RCSB PDB | REST | — | Generous. Search endpoint + per-entry fetch |
| STRING | REST | — | ~2 req/s. `required_score=700` for high confidence. `Semaphore(2)` |
| BioGRID | REST | `BIOGRID_ACCESS_KEY` | Moderate. Free registration at webservice.thebiogrid.org. Tool no-ops if key absent |
| DGIdb 5.0 | GraphQL | — | Generous. Dedupes drugs, surfaces approval status |
| ClinicalTrials.gov v2 | REST | — | Generous. Cloudflare-fronted; routed through `curl_cffi` Chrome 124 TLS impersonation in-client to handle WSL2/Cloudflare fingerprint blocks |
| Reactome | REST | — | ~5 req/s. Two-step token API. `Semaphore(3)` |
| SAbDab (OPIG) | REST / TSV | — | TSV download of full summary on first use; 7-day disk cache. RCSB FASTA + AbNum called for top-N structure CDR annotation |
| IEDB | PostgREST | — | Generous. B-cell epitope search by antigen description |
| IEDB NextGen Tools | REST async | — | Shared academic compute. Async ticket pattern (`/pipeline` POST → poll `/results/{id}`). 60s client-side timeout. NetMHCpan 4.1 EL recommended |
| gnomAD v4 | GraphQL | — | Generous. Constraint endpoint is fast; `gene.variants` payload is larger (session-cached per gene) |
| InterPro | REST | — | Generous. Keyed by UniProt accession |
| MaveDB | REST | — | Generous. Score-set search (JSON) + per-variant scores CSV. CSV cached per URN to avoid re-downloading for each variant |
| MyVariant.info | REST | — | Generous. Single HGVS genomic key returns ClinVar + dbNSFP (AlphaMissense, REVEL, CADD, SIFT, PolyPhen) + gnomAD AF |
| Ensembl REST + VEP | REST | — | 5 req/s (`Semaphore(5)`). HGNC → ENSG resolution, region overlap, VEP HGVS consequences (canonical transcript by default; `include_all_transcripts` opt-in) |
| GTEx | REST | — | Generous. Median TPM per tissue. Resolves HGNC → GENCODE via the shared Ensembl client |
| Human Protein Atlas | REST (bulk JSON) | — | Bulk-download JSON parsed for tissue-specificity category, subcellular localization, and pathology prognostics. 7-day disk cache |
| OpenFDA | REST | optional `OPENFDA_API_KEY` | Free tier: 240 req/min, 1000 req/day; key lifts the limit. FAERS adverse events, structured drug labels (boxed warnings), recall enforcement reports. `Semaphore(2)`. 7-day disk cache |
| NCBI E-utils | REST | — | 3 req/s. Set `NCBI_EMAIL` env var per NCBI ToS for production |

---

## Testing

```bash
uv run pytest tests/ -v          # 192 unit + integration tests
uv run pytest tests/ --cov=genesis_bio_mcp
```

Test layout:

- `tests/test_biochem.py` — ExPASy ProtParam reference-value tests for MW, pI, GRAVY, extinction, liability scanner (hen egg lysozyme, bovine RNase A)
- `tests/test_variant_parser.py` — protein-change parser (`R175H`, `p.Arg175His`, three-letter)
- `tests/test_clients.py` — per-client mocked HTTP tests (happy path, 404/empty, network error) using `respx` for httpx and `monkeypatch` for `curl_cffi`
- `tests/test_workflow_agent.py` — ToolSpec registry, agent loop, format_registry_docs
- `tests/test_e2e.py` — end-to-end integration tests with mocked full-pipeline runs

See [CONTRIBUTING.md](../CONTRIBUTING.md) for the full 7-step workflow to
add a new data source.
