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

        OT's ``search`` GraphQL is autocomplete-style and brittle to natural-
        language disease strings. ``"non-small-cell lung cancer (NSCLC)"``
        returns zero hits; even the cleaned ``"non-small-cell lung cancer"``
        returns zero. Only the bare abbreviation ``"NSCLC"`` resolves to
        ``EFO_0003060``. Failing silently here produces a confident-looking
        but wrong ``LOW PRIORITY`` call on real targets, so we explicitly
        retry with a small set of normalized variants and log which one won.
        """
        for candidate in _normalize_indication_variants(name):
            data = await self._graphql(_DISEASE_SEARCH_QUERY, {"name": candidate})
            if data is None:
                continue
            hits = data.get("data", {}).get("search", {}).get("hits", [])
            for hit in hits:
                if hit.get("entity") == "disease":
                    if candidate != name:
                        logger.info(
                            "Open Targets disease resolution: '%s' returned no hits, "
                            "fallback variant '%s' resolved to %s (%s)",
                            name,
                            candidate,
                            hit.get("name"),
                            hit["id"],
                        )
                    return hit["id"], hit.get("name")
        return None, None

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
