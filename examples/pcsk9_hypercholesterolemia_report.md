# Target Assessment: PCSK9 | hypercholesterolemia

**Priority Score: 9.1/10 — HIGH PRIORITY**

## Evidence Summary
PCSK9 shows strong Open Targets association with hypercholesterolemia (score: 0.83, n=4 evidence items). Open Targets reports strong known-drug evidence for PCSK9 (score: 0.98), suggesting existing approved or clinical-stage therapeutics — likely biologics if small-molecule data is sparse. DepMap CRISPR data show dependency in 5% of cancer lines. GWAS Catalog links 4 variants near PCSK9 to 'hypercholesterolemia'-related traits (strongest p=4.00e-20). ChEMBL reports 76 compounds with potency data against PCSK9; best IC50 ≈ 0.2 nM (clinical-grade, pChEMBL=9.8).

## Scoring Breakdown
| Source | Contribution | Max |
|---|---|---|
| Open Targets association | 2.48 | 3.0 |
| Cancer dependency | 0.1 (5% dependent) | 2.0 |
| GWAS evidence | 0.8 | 2.0 |
| Chemical matter | 1.5 (ChEMBL pChEMBL=9.8) | 1.5 |
| Protein annotation | 1.50 | 1.5 |

## Data Sources
- **UniProt:** ✓
- **Open Targets:** ✓
- **DepMap:** ✓
- **GWAS Catalog:** ✓
- **ChEMBL:** ✓
- **PubChem:** ✗ (no data)

**Data gaps:** pubchem

## Confidence Assessment
**Data coverage:** 83% of core sources returned data
**Score range:** 8.3–9.8/10 (uncertainty from 17% missing sources)

## API Latency
| API | Latency (s) |
|---|---|
| gwas | 15.02 ← slowest |
| chembl | 2.84 |
| pubchem | 1.94 |
| open_targets | 0.73 |
| depmap | 0.53 |
| uniprot | 0.18 |

---
_Resolved: PCSK9 | NCBI Gene: 255738 | UniProt: Q8NBP7_