# Target Assessment: CD274 | melanoma

**Priority Score: 4.1/10 — MEDIUM PRIORITY**

## Evidence Summary
CD274 shows modest Open Targets association with melanoma (score: 0.40, n=2 evidence items). Open Targets reports strong known-drug evidence for CD274 (score: 0.61), suggesting existing approved or clinical-stage therapeutics — likely biologics if small-molecule data is sparse. DepMap CRISPR data show dependency in 0% of cancer lines, highest in cutaneous melanoma, lymphoid neoplasm, multiple myeloma [DepMap score boosted 1.2× — indication matches top lineage]. ChEMBL reports 81 compounds with potency data against CD274; best IC50 ≈ 0.0 nM (clinical-grade, pChEMBL=10.7).

## Scoring Breakdown
| Source | Contribution | Max |
|---|---|---|
| Open Targets association | 1.2 | 3.0 |
| Cancer dependency | 0.0 (0% dependent (lineage match, 1.2×)) | 2.0 |
| GWAS evidence | 0.0 | 2.0 |
| Chemical matter | 1.5 (ChEMBL pChEMBL=10.7) | 1.5 |
| Protein annotation | 0.50 | 1.5 |

## Data Sources
- **UniProt:** ✓
- **Open Targets:** ✓
- **DepMap:** ✓
- **GWAS Catalog:** ✗ (no data)
- **ChEMBL:** ✓
- **PubChem:** ✗ (no data)

**Data gaps:** gwas, pubchem

## Confidence Assessment
**Data coverage:** 67% of core sources returned data
**Score range:** 3.4–4.8/10 (uncertainty from 33% missing sources)

## API Latency
| API | Latency (s) |
|---|---|
| pubchem | 1.70 ← slowest |
| gwas | 1.53 |
| chembl | 0.87 |
| open_targets | 0.52 |
| depmap | 0.36 |
| uniprot | 0.17 |

---
_Resolved: CD274 | NCBI Gene: 29126 | UniProt: Q9NZQ7_