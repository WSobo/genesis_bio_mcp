# Benchmark

End-to-end `prioritize_target` runs against live APIs — 12 representative
targets spanning oncology, metabolic disease, autoimmune, and pain biology
(single session, warm DepMap/EFO cache). This page tracks **coverage and
latency**; for how evidence is rated see
[architecture.md#evidence-profile-model--prioritize_target](architecture.md#evidence-profile-model--prioritize_target).

> **v0.7.0:** `prioritize_target` no longer emits a composite 0–10 score or
> High/Medium/Low tier. It returns a per-axis evidence profile instead, so the
> old `Score`/`Tier` columns have been retired from this table.

| Gene input | Resolved | Indication | GWAS | Time |
|---|---|---|---|---|
| BRAF | BRAF | melanoma | — somatic driver | 4.3s |
| EGFR | EGFR | non-small cell lung carcinoma | — | 2.2s |
| KRAS | KRAS | pancreatic cancer | — | 2.2s |
| **HER2** | **ERBB2** | breast cancer | — | 2.0s |
| PCSK9 | PCSK9 | hypercholesterolemia | ✓ 4 hits, p=4×10⁻²⁰ | 2.3s |
| FTO | FTO | obesity | timeout† | 14.4s |
| TNF | TNF | rheumatoid arthritis | ✓ 1 hit, p=9×10⁻²⁵ | 9.1s |
| PTGS2 | PTGS2 | inflammation | — | 7.2s |
| TP53 | TP53 | squamous cell carcinoma | — | 3.8s |
| CD274 | CD274 | melanoma | — | 6.3s |
| **p53** | **TP53** | lung cancer | — session cache | 3.3s |
| **COX2** | **PTGS2** | pain | — session cache | **1.9s** |

**Bold** gene inputs indicate alias resolution. GWAS gaps for BRAF / EGFR
/ KRAS are biologically expected — somatic cancer drivers are not GWAS
loci, so the genetics axis is rated from OT's `genetic_association` rather
than from absent GWAS hits. COX2/pain at 1.9s demonstrates the session gene
cache: PTGS2 associations were already fetched for the inflammation query
and are reused instantly.

†FTO has extensive GWAS signal (strong obesity associations with p < 10⁻¹⁰⁰)
but the Catalog API is slow for high-association genes; the 24h disk
cache rescues repeat queries.

---

## Example output — `prioritize_target("BRAF", "melanoma")`

### Standard mode — evidence profile

```
## Evidence Profile

| Axis                  | Signal       | Evidence |
|-----------------------|--------------|----------|
| Disease association   | ●●● Strong   | OT overall 0.82, n=5 evidence items |
| Clinical validation   | ●●● Strong   | known-drug 0.98 — approved / clinical-stage therapeutics exist |
| Genetic evidence      | ○○○ None     | OT genetic 0.00; no GWAS hits (somatic driver) |
| Cancer dependency     | ●●● Strong   | 9% of lines dependent; top: skin; matches indication |
| Chemical tractability | ●●● Strong   | best functional pChEMBL 9.5 (clinical-grade), 68 ChEMBL actives |
| Annotation quality    | ●●● Strong   | SwissProt-reviewed; 1 known variant |
```

Backing data the profile is read from:

```
→ disease_association.overall_score:          0.82
→ disease_association.somatic_mutation_score: 0.80   # BRAF V600E is the canonical somatic driver
→ disease_association.known_drug_score:       0.98   # vemurafenib, dabrafenib, encorafenib
→ cancer_dependency.fraction_dependent_lines: 0.09   # 9% — lineage-selective (melanoma/skin)
→ chembl_compounds.best_pchembl_functional:   9.5    # clinical-grade functional potency
→ data_gaps: ["gwas"]                                # expected — somatic, not germline
→ data_coverage_pct: 83.3                            # 5 of 6 core sources returned data
```

### Extended mode (`extended=True`)

Passes the same target through all four lab-loop tools in one parallel
gather and adds a **Tissue specificity** axis (HPA):

```
→ protein_structure.alphafold_plddt:    92.1       # high confidence (≥90)
→ protein_structure.best_resolution:    1.7 Å
→ protein_structure.has_ligand_bound:   true       # inhibitor co-crystal available
→ protein_interactome.top_partners:    MAP2K1 (0.999), MAP2K2 (0.998), RAF1 (0.963)
→ drug_history.approved_drug_count:    4
→ drug_history.trial_counts_by_phase:  {"Phase 1": 12, "Phase 2": 8, "Phase 3": 3}
→ pathway_context.top_pathway:         "MAPK1/MAPK2 Cascade" (p=2.3e-15)
```
