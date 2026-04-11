"""GWAS trait synonym mapping — fallback for EFO-backed matching.

Primary matching uses EFO ontology terms resolved via OLS4 (see efo_resolver.py).
This dict is the fallback when OLS4 is unavailable or returns nothing useful.

Each key is a canonical indication string. Values are GWAS Catalog free-text
label substrings known to cover that disease. Add entries here only when OLS4
consistently fails to resolve an indication that has real GWAS signal.

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

# Avoid circular import at type-check time
from typing import TYPE_CHECKING

from genesis_bio_mcp.models import GwasHit

if TYPE_CHECKING:
    from genesis_bio_mcp.config.efo_resolver import EFOTerm

# ---------------------------------------------------------------------------
# Fallback: canonical indication → GWAS Catalog trait label substrings
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
# Matching helpers
# ---------------------------------------------------------------------------


def _normalize(text: str) -> str:
    """Normalize Unicode to ASCII-comparable form."""
    return unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii").lower()


def filter_by_trait(
    hits: list[GwasHit],
    trait: str,
    efo_terms: list[EFOTerm] | None = None,
) -> list[GwasHit]:
    """Return hits whose trait matches the query indication.

    Matching strategy (in priority order):
    1. EFO URI exact match — precise, ontology-backed (gene-ID path only).
    2. String match against EFO labels + synonyms from OLS4 resolution.
    3. Fallback: string match against hardcoded TRAIT_SYNONYMS dict.

    The fallback ensures arbitrary indications still work even when OLS4 is
    unavailable, at the cost of requiring manual synonym maintenance.
    """
    trait_norm = _normalize(trait)

    # Build EFO URI set for precise matching (populated when efo_terms is provided)
    efo_uris: set[str] = set()
    match_strings: set[str] = {trait_norm}

    if efo_terms:
        for t in efo_terms:
            if t.uri:
                efo_uris.add(t.uri)
            match_strings.add(_normalize(t.label))
            match_strings.update(_normalize(s) for s in t.synonyms)
    else:
        # No EFO data — fall back to hardcoded synonyms
        match_strings.update(_normalize(s) for s in TRAIT_SYNONYMS.get(trait_norm, []))

    def _matches(hit: GwasHit) -> bool:
        # Precise: EFO URI match (gene-ID path associations embed efoTraits[].uri)
        if hit.efo_uri and hit.efo_uri in efo_uris:
            return True
        # Broader: label/synonym substring match (covers SNP path free-text labels)
        ht = _normalize(hit.trait)
        return any(term in ht for term in match_strings)

    return [h for h in hits if _matches(h)]
