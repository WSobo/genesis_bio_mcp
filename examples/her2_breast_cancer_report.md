# Target Assessment: ERBB2 | breast cancer

**Priority Score: 7.6/10 — HIGH PRIORITY**

## Evidence Summary
ERBB2 shows strong Open Targets association with breast cancer (score: 0.63, n=2 evidence items). Open Targets reports strong known-drug evidence for ERBB2 (score: 0.99), suggesting existing approved or clinical-stage therapeutics — likely biologics if small-molecule data is sparse. DepMap CRISPR data show dependency in 20% of cancer lines, highest in urinary bladder cancer, breast adenocarcinoma, urinary bladder carcinoma. ChEMBL reports 74 compounds with potency data against ERBB2; best IC50 ≈ 5.0 nM (lead-quality, pChEMBL=8.3).

## Scoring Breakdown
| Source | Contribution | Max |
|---|---|---|
| Open Targets association | 1.9 | 3.0 |
| Cancer dependency | 0.4 (20% dependent) | 2.0 |
| GWAS evidence | 0.0 | 2.0 |
| Chemical matter | 1.0 (ChEMBL pChEMBL=8.3) | 1.5 |
| Protein annotation | 1.50 | 1.5 |

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
**Score range:** 6.4–8.9/10 (uncertainty from 33% missing sources)

## API Latency
| API | Latency (s) |
|---|---|
| pubchem | 2.37 ← slowest |
| gwas | 1.72 |
| chembl | 1.50 |
| open_targets | 0.58 |
| depmap | 0.47 |
| uniprot | 0.24 |

---
_Resolved: ERBB2 | NCBI Gene: 2064 | UniProt: P04626_