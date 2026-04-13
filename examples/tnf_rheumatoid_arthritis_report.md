# Target Assessment: TNF | rheumatoid arthritis

**Priority Score: 6.5/10 — MEDIUM PRIORITY**

## Evidence Summary
TNF shows strong Open Targets association with rheumatoid arthritis (score: 0.64, n=3 evidence items). Open Targets reports strong known-drug evidence for TNF (score: 1.00), suggesting existing approved or clinical-stage therapeutics — likely biologics if small-molecule data is sparse. DepMap CRISPR data show dependency in 0% of cancer lines. ChEMBL reports 1 compounds with potency data against TNF; best IC50 ≈ 34.7 µM (hit-quality, pChEMBL=4.5).

## Scoring Breakdown
| Source | Contribution | Max |
|---|---|---|
| Open Targets association | 1.93 | 3.0 |
| Cancer dependency | 0.01 (0% dependent) | 2.0 |
| GWAS evidence | 0.0 | 2.0 |
| Chemical matter | 0.25 (ChEMBL pChEMBL=4.5) | 1.5 |
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
**Score range:** 6.0–7.0/10 (uncertainty from 17% missing sources)

## API Latency
| API | Latency (s) |
|---|---|
| gwas | 15.01 ← slowest |
| pubchem | 4.22 |
| chembl | 1.58 |
| open_targets | 0.59 |
| depmap | 0.58 |
| uniprot | 0.17 |

---
_Resolved: TNF | NCBI Gene: 7124 | UniProt: P01375_