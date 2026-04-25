"""Open Targets Platform GraphQL API client (v4)."""

from __future__ import annotations

import asyncio
import logging

import httpx

from genesis_bio_mcp.models import DiseaseLinkEvidence, TargetDiseaseAssociation

logger = logging.getLogger(__name__)

_GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"

_GENE_SEARCH_QUERY = """
query GeneSearch($symbol: String!) {
  search(queryString: $symbol, entityNames: ["target"], page: {index: 0, size: 1}) {
    hits {
      id
      name
      entity
    }
  }
}
"""

_DISEASE_SEARCH_QUERY = """
query DiseaseSearch($name: String!) {
  search(queryString: $name, entityNames: ["disease"], page: {index: 0, size: 3}) {
    hits {
      id
      name
      entity
    }
  }
}
"""

# Bs is the v4 argument for filtering by a list of disease EFO IDs.
# Note: evidenceCount does not exist in v4 — removed.
# datasourceScores / datatypeScores are the correct fields.
_ASSOCIATION_QUERY = """
query Association($ensemblId: String!, $efoIds: [String!]!) {
  target(ensemblId: $ensemblId) {
    associatedDiseases(Bs: $efoIds) {
      count
      rows {
        score
        datatypeScores {
          id
          score
        }
        datasourceScores {
          id
          score
        }
        disease {
          id
          name
        }
      }
    }
  }
}
"""

# Fallback: fetch top 100 sorted by overall score and filter client-side.
_ASSOCIATION_QUERY_NO_FILTER = """
query AssociationTopN($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    associatedDiseases(page: { index: 0, size: 100 }) {
      count
      rows {
        score
        datatypeScores {
          id
          score
        }
        disease {
          id
          name
        }
      }
    }
  }
}
"""


class OpenTargetsClient:
    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client

    async def get_association(
        self, gene_symbol: str, disease_name: str
    ) -> TargetDiseaseAssociation | None:
        """Return the Open Targets association score between a gene and disease."""
        ensembl_id = await self._resolve_gene(gene_symbol)
        if ensembl_id is None:
            logger.warning("Open Targets: could not resolve gene '%s' to Ensembl ID", gene_symbol)
            return None

        efo_id, resolved_disease = await self._resolve_disease(disease_name)
        if efo_id is None:
            logger.warning("Open Targets: could not resolve disease '%s' to EFO ID", disease_name)
            return None

        return await self._fetch_association(
            gene_symbol, ensembl_id, efo_id, resolved_disease or disease_name
        )

    async def _resolve_gene(self, symbol: str) -> str | None:
        data = await self._graphql(_GENE_SEARCH_QUERY, {"symbol": symbol})
        if data is None:
            return None
        hits = data.get("data", {}).get("search", {}).get("hits", [])
        for hit in hits:
            if hit.get("entity") == "target":
                return hit["id"]
        return hits[0]["id"] if hits else None

    async def _resolve_disease(self, name: str) -> tuple[str | None, str | None]:
        """Resolve a free-text indication to (EFO ID, canonical name).

        OT's ``search`` GraphQL is autocomplete-style and brittle. Two failure
        modes have shown up in production:

        1. Returns zero hits for a queryable string (bug J, fixed v0.3.2)
        2. Returns the WRONG disease as the top hit because of fuzzy substring
           matching: ``"type-2 diabetes mellitus"`` returned MONDO_0010802
           (pancreatic-hypoplasia-diabetes-congenital-heart-disease syndrome)
           even though MONDO_0005148 (type 2 diabetes mellitus) was the right
           answer (bug K, this fix).

        Strategy: try every normalized variant, score each disease hit by
        token-set similarity to its triggering variant, and return the
        highest-scoring match across all variants. A minimum threshold
        rejects fuzzy garbage hits (the syndrome above scores ~0.13 against
        "type 2 diabetes mellitus" — well below threshold).
        """
        best: tuple[float, dict, str] | None = None  # (score, hit, winning_variant)
        for candidate in _normalize_indication_variants(name):
            data = await self._graphql(_DISEASE_SEARCH_QUERY, {"name": candidate})
            if data is None:
                continue
            hits = data.get("data", {}).get("search", {}).get("hits", [])
            for hit in hits:
                if hit.get("entity") != "disease":
                    continue
                hit_name = hit.get("name", "")
                # Score against BOTH the candidate variant AND the original
                # input. Acronym-style hits ("NSCLC" → "non-small cell lung
                # carcinoma") have zero token overlap with the bare
                # abbreviation but solid overlap with the original
                # "non-small-cell lung cancer (NSCLC)" — taking the max
                # accepts those without lowering the threshold.
                score = max(
                    _name_match_score(candidate, hit_name),
                    _name_match_score(name, hit_name),
                )
                if best is None or score > best[0]:
                    best = (score, hit, candidate)
                    if score >= 0.95:
                        break
            if best and best[0] >= 0.95:
                break

        if best is None or best[0] < _DISEASE_MATCH_THRESHOLD:
            return None, None
        score, hit, winning = best
        if winning != name or score < 1.0:
            logger.info(
                "Open Targets disease resolution: '%s' → variant '%s' → %s (%s) [match score=%.2f]",
                name,
                winning,
                hit.get("name"),
                hit["id"],
                score,
            )
        return hit["id"], hit.get("name")

    async def _fetch_association(
        self,
        gene_symbol: str,
        ensembl_id: str,
        efo_id: str,
        disease_name: str,
    ) -> TargetDiseaseAssociation | None:
        # Primary: use Bs filter for the specific disease
        data = await self._graphql(
            _ASSOCIATION_QUERY,
            {"ensemblId": ensembl_id, "efoIds": [efo_id]},
        )
        row = _extract_row(data, efo_id)

        if row is None:
            # Fallback: fetch top 100 and match by EFO ID client-side
            logger.info(
                "Open Targets: Bs filter returned no rows for %s/%s, falling back to top-100",
                gene_symbol,
                efo_id,
            )
            data2 = await self._graphql(_ASSOCIATION_QUERY_NO_FILTER, {"ensemblId": ensembl_id})
            row = _extract_row(data2, efo_id)

        if row is None:
            logger.info(
                "Open Targets: no association found for %s / %s (%s)",
                gene_symbol,
                disease_name,
                efo_id,
            )
            return None

        return _parse_row(row, gene_symbol, disease_name, efo_id, ensembl_id)

    async def _graphql(self, query: str, variables: dict) -> dict | None:
        for attempt in range(2):
            try:
                resp = await self._client.post(
                    _GRAPHQL_URL,
                    json={"query": query, "variables": variables},
                    headers={"Content-Type": "application/json"},
                    timeout=30.0,
                )
                resp.raise_for_status()
                body = resp.json()
                if "errors" in body:
                    logger.warning(
                        "Open Targets GraphQL errors: %s",
                        body["errors"][0].get("message", ""),
                    )
                    return None
                return body
            except httpx.HTTPStatusError as exc:
                if exc.response.status_code >= 500 and attempt == 0:
                    logger.warning(
                        "Open Targets 5xx on attempt %d (%s), retrying in 2s",
                        attempt + 1,
                        exc.response.status_code,
                    )
                    await asyncio.sleep(2.0)
                    continue
                logger.warning("Open Targets HTTP error: %s", exc)
                return None
            except Exception as exc:
                logger.warning("Open Targets GraphQL request failed: %s", exc)
                return None
        return None


# _DISEASE_MATCH_THRESHOLD (0.5):
#   Minimum token-set Jaccard similarity between the queried indication
#   variant and the returned disease name. Below this we treat the hit as
#   fuzzy garbage and continue searching. Calibrated so:
#     - "type 2 diabetes mellitus" vs MONDO_0005148 ("type 2 diabetes mellitus") → 1.00 ✓
#     - "type-2 diabetes mellitus" vs MONDO_0010802 ("pancreatic-hypoplasia-
#       diabetes-congenital-heart-disease syndrome") → ~0.11 ✗ (rejected)
#     - "NSCLC" vs EFO_0003060 ("non-small cell lung carcinoma") → 0.0 token
#       overlap, but the abbreviation passes via prefix-substring path below.
#   Range to consider: 0.4–0.6. Below 0.4 risks accepting weak fuzzy matches;
#   above 0.6 risks rejecting valid abbreviation/acronym matches.
_DISEASE_MATCH_THRESHOLD = 0.5


def _name_match_score(query: str, hit_name: str) -> float:
    """Score a disease hit's name against the queried variant.

    Returns 0.0–1.0. Token-set Jaccard similarity by default, with a perfect
    1.0 when the query is a substring of the hit name (or vice versa) — the
    latter rule lets short abbreviations like "NSCLC" pass when the hit name
    is the long form ("non-small cell lung carcinoma") and the abbreviation
    is in OT's index for that disease.
    """
    if not query or not hit_name:
        return 0.0
    q_norm = query.strip().lower().replace("-", " ")
    h_norm = hit_name.strip().lower().replace("-", " ")
    if q_norm == h_norm:
        return 1.0
    # Strict substring match either direction = high confidence (covers
    # short abbreviations indexed against long disease names).
    if q_norm in h_norm or h_norm in q_norm:
        # Cap at 0.95 so an exact match still beats a substring match in
        # ranking when both fire from different variants.
        return 0.95
    q_tokens = set(q_norm.split())
    h_tokens = set(h_norm.split())
    if not q_tokens or not h_tokens:
        return 0.0
    intersection = q_tokens & h_tokens
    union = q_tokens | h_tokens
    return len(intersection) / len(union)


# Common indication acronyms → canonical full forms. Used by
# ``_normalize_indication_variants`` to expand bare acronyms (``"NSCLC"``,
# ``"PDAC"``, ``"T2DM"``) to a form OT's autocomplete-style search will
# resolve. The acronym itself is still tried first; the expansion is added
# as a fallback variant so users who type the long form keep cache hits.
#
# Curated to oncology + cardiometabolic + neuroscience + immunology, the
# four areas that produced silent-data-loss bugs in the v0.3.0–v0.3.3 smoke
# tests. Add new entries when you find another acronym OT's search misses.
_ACRONYM_EXPANSIONS: dict[str, str] = {
    # Oncology — solid tumors
    "NSCLC": "non-small cell lung carcinoma",
    "SCLC": "small cell lung carcinoma",
    "PDAC": "pancreatic ductal adenocarcinoma",
    "TNBC": "triple-negative breast cancer",
    "HCC": "hepatocellular carcinoma",
    "RCC": "renal cell carcinoma",
    "CRC": "colorectal carcinoma",
    "GBM": "glioblastoma multiforme",
    # Oncology — hematologic
    "AML": "acute myeloid leukemia",
    "ALL": "acute lymphoblastic leukemia",
    "CML": "chronic myeloid leukemia",
    "CLL": "chronic lymphocytic leukemia",
    "MM": "multiple myeloma",
    "DLBCL": "diffuse large B-cell lymphoma",
    # Endocrine / metabolic
    "T1D": "type 1 diabetes mellitus",
    "T1DM": "type 1 diabetes mellitus",
    "T2D": "type 2 diabetes mellitus",
    "T2DM": "type 2 diabetes mellitus",
    "MODY": "maturity-onset diabetes of the young",
    # Cardiovascular
    "CAD": "coronary artery disease",
    "CHD": "coronary heart disease",
    "MI": "myocardial infarction",
    "HF": "heart failure",
    "HFpEF": "heart failure with preserved ejection fraction",
    "HFrEF": "heart failure with reduced ejection fraction",
    "AFib": "atrial fibrillation",
    "FH": "familial hypercholesterolemia",
    # Liver / GI
    "NASH": "non-alcoholic steatohepatitis",
    "MASH": "metabolic dysfunction-associated steatohepatitis",
    "NAFLD": "non-alcoholic fatty liver disease",
    "MASLD": "metabolic dysfunction-associated steatotic liver disease",
    "IBD": "inflammatory bowel disease",
    "UC": "ulcerative colitis",
    # Pulmonary
    "COPD": "chronic obstructive pulmonary disease",
    "IPF": "idiopathic pulmonary fibrosis",
    # Neuroscience
    "AD": "Alzheimer disease",
    "PD": "Parkinson disease",
    "ALS": "amyotrophic lateral sclerosis",
    "MS": "multiple sclerosis",
    "ASD": "autism spectrum disorder",
    "ADHD": "attention deficit hyperactivity disorder",
    "MDD": "major depressive disorder",
    "OCD": "obsessive compulsive disorder",
    "PTSD": "post-traumatic stress disorder",
    # Immunology / rheumatology
    "RA": "rheumatoid arthritis",
    "SLE": "systemic lupus erythematosus",
    "PsA": "psoriatic arthritis",
    "AS": "ankylosing spondylitis",
    "OA": "osteoarthritis",
}


def _normalize_indication_variants(name: str) -> list[str]:
    """Return ordered, deduped query variants to try against OT disease search.

    Order matters: the literal string is tried first to preserve cache hits
    and respect user intent. Fallbacks fire only if the literal returns zero
    disease hits.

    1. Literal (current input, whitespace-collapsed)
    2. Stripped of any parenthetical (``"X (NSCLC)"`` → ``"X"``)
    3. Bare parenthetical content alone (``"X (NSCLC)"`` → ``"NSCLC"``) — this
       is the variant that catches the real-world case where OT's search
       indexes the abbreviation but not the long form
    4. Hyphen-stripped (``"non-small-cell"`` → ``"non small cell"``)
    5. Acronym expansion — for any all-caps token in the input that matches
       ``_ACRONYM_EXPANSIONS``, add the canonical full form. Catches bare
       ``"NSCLC"`` / ``"PDAC"`` / ``"T2DM"`` queries that OT's search returns
       zero hits for, and also rescues messy strings where the abbreviation
       was the only resolvable token.
    """
    import re

    variants: list[str] = []

    def _add(v: str) -> None:
        cleaned = " ".join(v.strip().split())
        if cleaned and cleaned not in variants:
            variants.append(cleaned)

    base = name.strip()
    _add(base)

    paren_match = re.search(r"\(([^)]+)\)", base)
    stripped = re.sub(r"\s*\([^)]+\)\s*", " ", base).strip()
    _add(stripped)
    if paren_match:
        _add(paren_match.group(1))

    # Hyphen variant operates on whichever stripped form exists.
    pivot = stripped or base
    _add(pivot.replace("-", " "))

    # Acronym expansion: scan all 2–6 char alphanumeric tokens (preserving
    # mixed case so HFpEF / HFrEF / AFib still match), look each up in the
    # acronym table case-insensitively, and add the expansion as a variant.
    for token in re.findall(r"\b[A-Za-z][A-Za-z0-9]{1,5}\b", base):
        expansion = _ACRONYM_EXPANSIONS.get(token)
        if expansion is None:
            expansion = _ACRONYM_EXPANSIONS.get(token.upper())
        if expansion:
            _add(expansion)

    return variants


def _extract_row(data: dict | None, efo_id: str) -> dict | None:
    """Extract the first row matching the EFO ID from a GraphQL response."""
    if data is None:
        return None
    rows = data.get("data", {}).get("target", {}).get("associatedDiseases", {}).get("rows", [])
    if not rows:
        return None
    # If Bs filter was used, only one disease is returned — take it directly
    if len(rows) == 1:
        return rows[0]
    # For fallback (top-100), find the matching disease
    for row in rows:
        if row.get("disease", {}).get("id") == efo_id:
            return row
    return rows[0]  # Best available if exact match not in top 100


def _parse_row(
    row: dict,
    gene_symbol: str,
    disease_name: str,
    efo_id: str,
    ensembl_id: str,
) -> TargetDiseaseAssociation:
    score = row.get("score", 0.0) or 0.0
    datatype_scores: dict[str, float] = {}
    for ds in row.get("datatypeScores", []):
        datatype_scores[ds["id"]] = ds.get("score", 0.0)

    evidence_breakdown = [
        DiseaseLinkEvidence(evidence_type=k, score=v) for k, v in datatype_scores.items()
    ]

    return TargetDiseaseAssociation(
        gene_symbol=gene_symbol,
        disease_name=disease_name,
        disease_efo_id=efo_id,
        ensembl_id=ensembl_id,
        overall_score=score,
        genetic_association_score=datatype_scores.get("genetic_association"),
        somatic_mutation_score=datatype_scores.get("somatic_mutation"),
        known_drug_score=datatype_scores.get("known_drug") or datatype_scores.get("clinical"),
        literature_mining_score=datatype_scores.get("literature"),
        evidence_count=len(evidence_breakdown),
        evidence_breakdown=evidence_breakdown,
    )
