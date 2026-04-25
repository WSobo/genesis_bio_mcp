# Tool catalog

genesis-bio-mcp exposes **27 MCP tools** across **23 biomedical data sources**
(plus `run_biology_workflow`, the 28th, which chains the others with Claude).

All tools return Markdown strings by default and accept
`response_format="json"` for programmatic integration. Single-gene tools
auto-resolve aliases via `resolve_gene` (`HER2` → `ERBB2`, `p53` → `TP53`,
`COX2` → `PTGS2`).

Jump to:
[Gene annotation](#gene-annotation) ·
[Disease evidence](#disease-evidence) ·
[Druggability](#druggability) ·
[Structure & interactions](#structure--interactions) ·
[Antibody & epitope](#antibody--epitope) ·
[Protein engineering](#protein-engineering) ·
[Expression](#expression) ·
[Pathways](#pathways) ·
[Synthesis & orchestration](#synthesis--orchestration)

---

## Gene annotation

| Tool | Data source | Purpose |
|---|---|---|
| `resolve_gene` | UniProt + NCBI | Resolve aliases → canonical HGNC symbol + NCBI Gene ID + UniProt accession |
| `get_protein_info` | UniProt Swiss-Prot | Function summary, subcellular location, disease variants, PDB cross-refs, disulfide bonds |
| `get_protein_sequence` | UniProt FASTA | Protein sequence (optional residue slice) + biochem features + liability-motif scan — see [Protein engineering](#protein-engineering) |

## Disease evidence

| Tool | Data source | Purpose |
|---|---|---|
| `get_target_disease_association` | Open Targets v4 GraphQL | Evidence-based association score (0–1) + per-datatype breakdown (genetic, somatic, known_drug, literature) |
| `get_cancer_dependency` | DepMap + Open Targets | CRISPR essentiality fraction, pan-essential flag, top dependent lineages. 7-day disk cache, OT lineage proxy fallback |
| `get_gwas_evidence` | GWAS Catalog + EFO (OLS4) | Genome-wide significant SNP associations with EFO-backed trait matching |

## Druggability

| Tool | Data source | Purpose |
|---|---|---|
| `get_compounds` | PubChem | Active small molecules with bioactivity hit counts (quick tractability check) |
| `get_chembl_compounds` | ChEMBL | Quantitative potency (IC50 / Ki / Kd) with pChEMBL values and confidence scores |

## Structure & interactions

| Tool | Data source | Purpose |
|---|---|---|
| `get_protein_structure` | AlphaFold + RCSB PDB | pLDDT confidence, experimental resolution, ligand-bound structures |
| `get_protein_interactome` | STRING | High-confidence (score ≥700) binding partners with evidence-channel breakdown |
| `get_biogrid_interactions` | BioGRID | Curated literature PPI with experimental method annotation. Requires `BIOGRID_ACCESS_KEY` env var |

## Antibody & epitope

| Tool | Data source | Purpose |
|---|---|---|
| `get_antibody_structures` | SAbDab + RCSB FASTA + AbNum | Curated antibody/nanobody crystal structures with CDR annotation (Chothia scheme) for top results |
| `get_epitope_data` | IEDB | Known B-cell epitope records with isotypes, linear sequences, and PDB structural evidence |
| `get_mhc_binding` | IEDB NextGen Tools (NetMHCpan 4.1) | T-cell epitope / MHC-I/II binding prediction — see [Protein engineering](#protein-engineering) |

## Protein engineering

Introduced in v0.2.0 — these are the primary tools for protein/antibody
engineering workflows. See [docs/protein-engineering.md](protein-engineering.md)
for end-to-end examples.

| Tool | Data source | Purpose |
|---|---|---|
| `get_protein_sequence` | UniProt FASTA + pure-Python biochem | FASTA + molecular weight, theoretical pI, GRAVY, net charge, ε₂₈₀ reduced/oxidized + liability-motif scan (deamidation NG/NS, isomerization DG/DS, N-glycosylation NXS/NXT, oxidation M/W, free cysteines via UniProt DISULFID annotation) |
| `get_variant_effects` | gnomAD + MyVariant.info + MaveDB + Ensembl VEP | Mutation-level pathogenicity: ClinVar submissions, AlphaMissense class + score, REVEL, CADD-Phred, SIFT, PolyPhen-2, gnomAD allele frequency, DMS fitness scores from MaveDB, plus VEP splice/UTR/regulatory consequences |
| `get_variant_constraints` | gnomAD v4 | Gene-level pLI, LOEUF, oe_lof, oe_mis, Z-scores — tolerance-to-mutation pre-filter |
| `get_variant_consequences` | Ensembl VEP | Splice / UTR / regulatory overlap and SIFT / PolyPhen for any HGVS or coordinate variant. Canonical-transcript by default; `include_all_transcripts` opt-in. Complements `get_variant_effects` (dbNSFP scores) with VEP's transcript-aware consequence calls |
| `get_domain_annotation` | InterPro | Domain boundaries with Pfam/SMART/CDD/GO term annotation |
| `get_dms_scores` | MaveDB | Available deep mutational scanning score-set metadata (URNs, variant counts, citations) |
| `get_mhc_binding` | IEDB NextGen Tools | MHC-I / II binding predictions for peptides or auto-windowed proteins against a configurable HLA panel |

## Expression

Tissue and protein expression evidence — added in v0.3.0 to support
indication-aware target prioritization (the `expression` scoring axis pulls
from `get_protein_atlas`).

| Tool | Data source | Purpose |
|---|---|---|
| `get_tissue_expression` | GTEx v8 | Bulk-RNA median TPM across ~54 human tissues. Resolves HGNC → GENCODE via the shared Ensembl client. Within `prioritize_target` extended mode, the report is restricted to therapeutic-area-relevant tissues via `config/indication_tissue_map.py` |
| `get_protein_atlas` | Human Protein Atlas | Tissue-specificity category (`Tissue enriched`, `Group enriched`, `Tissue enhanced`, `Low tissue specificity`, `Not detected`), subcellular localization, and pathology prognostics (cancer-survival hazard ratios). Powers the v0.3.0 expression scoring axis |

## Pathways

| Tool | Data source | Purpose |
|---|---|---|
| `get_pathway_context` | Reactome | Pathway enrichment for a gene (FDR-ranked pathway membership) |
| `get_pathway_members` | Reactome | Enumerate all HGNC gene symbols in a named pathway (display name or stable ID) |

## Synthesis & orchestration

| Tool | Data source | Purpose |
|---|---|---|
| `get_drug_history` | DGIdb + ClinicalTrials.gov + OpenFDA | Known drug interactions with approval status + trial counts by phase. Top-5 approved drugs are enriched with OpenFDA post-market safety signals (FAERS adverse-event counts, boxed warnings, recalls) carrying a permanent regulatory disclaimer |
| `prioritize_target` | All of the above in parallel | Full evidence synthesis → composite priority score (0–10) with confidence interval. See [docs/architecture.md](architecture.md) for the scoring model |
| `compare_targets` | Calls `prioritize_target` ×N | Side-by-side comparison of 2–5 targets for one indication |
| `run_biology_workflow` | Claude + ToolSpec registry | Free-text question → inner Claude agent (`claude-sonnet-4-6`) dynamically selects and chains tools. Requires `ANTHROPIC_API_KEY` |

---

## Common input patterns

Every single-gene tool accepts `gene_symbol` (HGNC preferred, aliases auto-resolved)
and `response_format` (`"markdown"` or `"json"`).

```python
# Minimum input
{"gene_symbol": "BRAF"}

# JSON output for pipeline integration
{"gene_symbol": "BRAF", "response_format": "json"}
```

Exceptions (tools with additional required fields):

| Tool | Extra required fields |
|---|---|
| `resolve_gene` | `gene_name` (instead of `gene_symbol`; accepts aliases) |
| `get_target_disease_association` | `disease_name` |
| `get_gwas_evidence` | `trait` |
| `get_protein_sequence` | optional `start`, `end` (1-indexed inclusive slice) |
| `get_variant_effects` | `mutation` (e.g. `R175H`, `p.Arg175His`) |
| `get_variant_consequences` | one of: `gene_symbol` + `mutation`; `hgvs_genomic` (e.g. `7:g.140753336A>T`); or `chrom` + `pos` + `ref` + `alt`. Optional `include_all_transcripts` |
| `get_antibody_structures` | `antigen_query` (protein name for best results), optional `max_results` |
| `get_epitope_data` | `antigen_query` |
| `get_mhc_binding` | `sequence`; optional `hla_alleles`, `mhc_class`, `peptide_lengths`, `method` |
| `get_pathway_members` | `pathway_name_or_id`, optional `max_genes` |
| `prioritize_target` / `compare_targets` | `indication` (+ `gene_symbols` list for compare); optional `extended=True` |
| `run_biology_workflow` | `question` (free text) |

---

## MCP Resource: `tool://registry`

Any MCP client can read the `tool://registry` resource to get a structured
Markdown catalogue of all 28 tools grouped by category, with the
embedding-searchable `use_when` guidance the workflow agent uses for
semantic retrieval. No tool call needed — useful for agent frameworks doing
capability auditing or embedding-based tool selection.
