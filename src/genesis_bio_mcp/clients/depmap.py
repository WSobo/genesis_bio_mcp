"""Cancer dependency client.

DepMap does not expose a stable public REST API (all portal endpoints return 404).
This client queries the Open Targets Platform for cancer-specific somatic mutation
evidence across cancer types, which reflects the same cancer dependency signal that
DepMap CRISPR screens capture (BRAF is essential in melanoma/lung/thyroid because
these cancers have activating BRAF mutations and are 'oncogene addicted').

The returned CancerDependency model is populated with:
  - Somatic mutation scores per cancer type (from Open Targets datasourceScores)
  - fraction_dependent_lines: fraction of cancer types with strong somatic evidence (>0.5)
  - pan_essential: False (CRISPR pan-essentiality is not inferrable from OT data)
  - data_source: clearly labelled as Open Targets somatic mutation evidence
"""

from __future__ import annotations

import logging
from typing import Optional

import httpx

from genesis_bio_mcp.models import CancerDependency, CellLineEssentiality

logger = logging.getLogger(__name__)

_GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"
_DEPENDENCY_THRESHOLD = 0.5   # OT somatic mutation score threshold (not CERES)

# Query: get all cancer-type associations sorted by somatic_mutation score
_CANCER_ASSOCIATIONS_QUERY = """
query CancerEvidence($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    associatedDiseases(orderByScore: "somatic_mutation", page: { index: 0, size: 50 }) {
      count
      rows {
        score
        datatypeScores { id score }
        disease {
          id
          name
          therapeuticAreas { name }
        }
      }
    }
  }
}
"""

# Gene symbol → Ensembl ID (reuse same pattern as open_targets.py)
_GENE_SEARCH_QUERY = """
query GeneSearch($symbol: String!) {
  search(queryString: $symbol, entityNames: ["target"], page: {index: 0, size: 1}) {
    hits { id name entity }
  }
}
"""


class DepMapClient:
    """Cancer dependency client backed by Open Targets somatic mutation evidence."""

    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client

    async def get_essentiality(self, gene_symbol: str) -> Optional[CancerDependency]:
        """Return cancer dependency evidence for a gene using Open Targets somatic data."""
        symbol = gene_symbol.strip().upper()

        ensembl_id = await self._resolve_gene(symbol)
        if ensembl_id is None:
            logger.warning("DepMap/OT: could not resolve gene '%s' to Ensembl ID", symbol)
            return None

        return await self._fetch_cancer_evidence(symbol, ensembl_id)

    async def _resolve_gene(self, symbol: str) -> Optional[str]:
        data = await self._graphql(_GENE_SEARCH_QUERY, {"symbol": symbol})
        if data is None:
            return None
        hits = data.get("data", {}).get("search", {}).get("hits", [])
        for hit in hits:
            if hit.get("entity") == "target":
                return hit["id"]
        return hits[0]["id"] if hits else None

    async def _fetch_cancer_evidence(
        self, gene_symbol: str, ensembl_id: str
    ) -> Optional[CancerDependency]:
        data = await self._graphql(_CANCER_ASSOCIATIONS_QUERY, {"ensemblId": ensembl_id})
        if data is None:
            return None

        rows = (
            data.get("data", {})
            .get("target", {})
            .get("associatedDiseases", {})
            .get("rows", [])
        )
        if not rows:
            return None

        # Filter to cancer-related diseases
        cancer_rows = [
            r for r in rows
            if _is_cancer(r)
        ]
        if not cancer_rows:
            cancer_rows = rows  # Fall back to all if no cancer filter match

        # Extract somatic mutation score per cancer type as proxy for dependency
        cell_lines: list[CellLineEssentiality] = []
        for row in cancer_rows:
            sm_score = next(
                (d["score"] for d in row.get("datatypeScores", []) if d["id"] == "somatic_mutation"),
                0.0,
            )
            if sm_score == 0.0:
                continue
            disease_name = row["disease"]["name"]
            areas = [a["name"] for a in row["disease"].get("therapeuticAreas", [])]
            lineage = next((a for a in areas if "cancer" in a.lower() or "tumor" in a.lower()), areas[0] if areas else "Oncology")
            # Map OT somatic score to CERES-like scale: invert (high OT score → negative CERES proxy)
            ceres_proxy = -(sm_score * 2.0)  # Range approx 0 to -2
            cell_lines.append(
                CellLineEssentiality(
                    cell_line=disease_name[:50],
                    lineage=lineage,
                    ceres_score=round(ceres_proxy, 3),
                    is_dependent=sm_score >= _DEPENDENCY_THRESHOLD,
                )
            )

        if not cell_lines:
            return None

        scores = [cl.ceres_score for cl in cell_lines]
        mean_score = sum(scores) / len(scores)
        n_dependent = sum(1 for cl in cell_lines if cl.is_dependent)
        fraction_dependent = n_dependent / len(cell_lines)

        # Top lineages (most dependent cancer types)
        lineage_groups: dict[str, list[float]] = {}
        for cl in cell_lines:
            lineage_groups.setdefault(cl.lineage, []).append(cl.ceres_score)
        lineage_means = {k: sum(v) / len(v) for k, v in lineage_groups.items()}
        top_lineages = sorted(lineage_means, key=lambda k: lineage_means[k])[:5]

        top_lines = sorted(cell_lines, key=lambda cl: cl.ceres_score)[:10]

        return CancerDependency(
            gene_symbol=gene_symbol,
            mean_ceres_score=round(mean_score, 4),
            fraction_dependent_lines=round(fraction_dependent, 4),
            pan_essential=False,  # Not determinable from OT data
            top_dependent_lineages=top_lineages,
            cell_lines=top_lines,
            data_source=(
                "Open Targets Platform v4 — somatic mutation evidence across cancer types "
                "(proxy for cancer dependency; CERES scores are scaled OT somatic mutation scores)"
            ),
        )

    async def _graphql(self, query: str, variables: dict) -> Optional[dict]:
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
                logger.warning("OT cancer query errors: %s", body["errors"][0].get("message", ""))
                return None
            return body
        except Exception as exc:
            logger.warning("DepMap/OT query failed: %s", exc)
            return None


def _is_cancer(row: dict) -> bool:
    areas = [a["name"].lower() for a in row["disease"].get("therapeuticAreas", [])]
    return any("cancer" in a or "tumor" in a or "neoplasm" in a or "oncology" in a for a in areas)
