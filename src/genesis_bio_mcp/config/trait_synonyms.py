"""GWAS trait synonym mapping for query–association matching.

Each key is the canonical indication string passed to the GWAS client.
Values are the trait label substrings that GWAS Catalog uses for that disease
(GWAS Catalog free-text trait labels are inconsistent; EFO ID-based matching
is the long-term upgrade path once EFO term lookup is wired in).

Maintainers: add synonyms here when a new indication returns a GWAS gap despite
known associations existing in the Catalog. Use GWAS Catalog's Trait search
(https://www.ebi.ac.uk/gwas/search) to find the exact label variants.

Design notes:
- Matching is case-insensitive substring: "waist" catches "waist circumference"
  and "waist-hip ratio" without needing them spelled out.
- Prefer specific substrings over short words: "adipose" yes, "fat" no
  (too many false positives in cancer/inflammation contexts).
- Both US and UK spellings should be listed where GWAS Catalog is inconsistent
  (e.g. "hypercholesterolemia" / "hypercholesterolaemia").
"""

from __future__ import annotations

import unicodedata

from genesis_bio_mcp.models import GwasHit

# ---------------------------------------------------------------------------
# Canonical indication → GWAS Catalog trait label substrings
# ---------------------------------------------------------------------------

TRAIT_SYNONYMS: dict[str, list[str]] = {
    # Lipid / PCSK9 / statin targets — GWAS Catalog uses both the clinical term
    # and the biomarker names (LDL, total cholesterol, lipid levels).
    "hypercholesterolemia": ["cholesterol", "ldl", "low-density lipoprotein", "lipid"],
    "hypercholesterolaemia": ["cholesterol", "ldl", "low-density lipoprotein", "lipid"],
    # Obesity / FTO — GWAS Catalog labels vary widely. "body mass index" and
    # "obesity" are most common; "adipose"/"fat mass"/"body weight" appear for
    # body composition GWAS. "waist" catches waist circumference and waist-hip ratio.
    # "fat" alone is excluded — too broad (cancer fat-mass confounders, etc.).
    "obesity": [
        "obesity",
        "body mass index",
        "bmi",
        "adiposity",
        "adipose",
        "overweight",
        "waist",
        "body weight",
        "body fat",
        "fat mass",
    ],
    # Autoimmune — "arthritis" alone is too broad (includes OA, psoriatic);
    # keeping the compound term as the primary.
    "rheumatoid arthritis": ["rheumatoid arthritis", "arthritis"],
    # Inflammation markers — CRP is the canonical biomarker for acute inflammation
    # GWAS; PTGS2/COX-2 studies often appear under CRP or "inflammatory marker".
    "inflammation": ["inflammation", "inflammatory", "c-reactive protein", "crp"],
    # Cardiovascular — covers CAD, MI, CHD under a single query.
    "cardiovascular disease": [
        "cardiovascular",
        "coronary artery",
        "myocardial infarction",
        "heart disease",
    ],
    # T2D — "glycated haemoglobin" / HbA1c studies are the predominant GWAS
    # signal for T2D susceptibility loci; including both spellings.
    "type 2 diabetes": [
        "type 2 diabetes",
        "t2d",
        "diabetes mellitus",
        "glycated haemoglobin",
        "hba1c",
    ],
    # Neurodegeneration — Alzheimer's GWAS labels vary; "dementia" catches
    # studies that don't specify subtype.
    "alzheimer disease": ["alzheimer", "dementia", "cognitive decline"],
    # Oncology — lung histology terms used by GWAS Catalog.
    "non-small cell lung carcinoma": [
        "lung cancer",
        "lung carcinoma",
        "non-small cell lung",
    ],
    "squamous cell carcinoma": ["squamous cell", "carcinoma"],
    # Pain — PTGS2/COX-2; "pain" is broad but matches the indication directly.
    "pain": ["pain"],
}


# ---------------------------------------------------------------------------
# Matching helper
# ---------------------------------------------------------------------------


def _normalize(text: str) -> str:
    """Normalize Unicode to ASCII-comparable form."""
    return unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii").lower()


def filter_by_trait(hits: list[GwasHit], trait: str) -> list[GwasHit]:
    """Return hits whose trait label matches the query indication.

    Falls back to direct substring match when the indication has no synonym entry,
    so arbitrary indications still work (just without synonym expansion).
    """
    trait_norm = _normalize(trait)
    synonyms = TRAIT_SYNONYMS.get(trait_norm, [])
    match_terms = [trait_norm] + [_normalize(s) for s in synonyms]

    def _matches(hit_trait: str) -> bool:
        ht = _normalize(hit_trait)
        return any(term in ht for term in match_terms)

    return [h for h in hits if _matches(h.trait)]
