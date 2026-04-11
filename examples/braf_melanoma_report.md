# Target Assessment: BRAF | melanoma

**Priority Score: 7.2/10 — HIGH PRIORITY**

## Evidence Summary
BRAF shows strong Open Targets association with melanoma (score: 0.82, n=5 evidence items). Open Targets reports strong known-drug evidence for BRAF (score: 0.98), suggesting existing approved or clinical-stage therapeutics — likely biologics if small-molecule data is sparse. DepMap CRISPR data show dependency in 9% of cancer lines, highest in differentiated thyroid carcinoma, glioblastoma multiforme, lung adenocarcinoma [DepMap score boosted 1.2× — indication matches top lineage]. ChEMBL reports 68 compounds with potency data against BRAF; best IC50 ≈ 0.3 nM (clinical-grade, pChEMBL=9.5).

## Scoring Breakdown
| Source | Contribution | Max |
|---|---|---|
| Open Targets association | 2.46 | 3.0 |
| Cancer dependency | 0.23 (9% dependent (lineage match, 1.2×)) | 2.0 |
| GWAS evidence | 0.0 | 2.0 |
| Chemical matter | 1.5 (ChEMBL pChEMBL=9.5) | 1.5 |
| Protein annotation | 1.50 | 1.5 |

## Data Sources
- **UniProt:** ✓
- **Open Targets:** ✓
- **DepMap:** ✓
- **GWAS Catalog:** ✗ (no data)
- **ChEMBL:** ✓
- **PubChem:** ✓

**Data gaps:** gwas

## Confidence Assessment
**Data coverage:** 83% of core sources returned data
**Score range:** 6.5–7.8/10 (uncertainty from 17% missing sources)

## API Latency
| API | Latency (s) |
|---|---|
| gwas | 3.22 ← slowest |
| pubchem | 2.19 |
| chembl | 2.14 |
| open_targets | 0.66 |
| depmap | 0.50 |
| uniprot | 0.19 |

---
_Resolved: BRAF | NCBI Gene: 673 | UniProt: P15056_