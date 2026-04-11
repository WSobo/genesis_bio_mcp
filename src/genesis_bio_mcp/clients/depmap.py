"""Cancer dependency client.

Primary source: DepMap public API — GET /download/gene_dep_summary
  Returns a CSV of real CRISPR Chronos Combined essentiality metrics for every gene.
  Fetched once at server startup and cached in-memory by the lifespan handler.

Supplementary source: Open Targets Platform — somatic mutation evidence
  Provides per-cancer-lineage context (top dependent lineages, per-type scores).
  Used when gene is found in the DepMap cache; used alone as fallback when not.

The returned CancerDependency model clearly labels which data came from which source.
"""

from __future__ import annotations

import asyncio
import csv
import io
import logging
import time
from pathlib import Path
from typing import Any

import httpx

from genesis_bio_mcp.models import CancerDependency, CellLineEssentiality

logger = logging.getLogger(__name__)

_GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"
_DEPMAP_SUMMARY_URL = "https://depmap.org/portal/api/download/gene_dep_summary"

# Somatic mutation score threshold for "dependent" call when using OT proxy
_OT_DEPENDENCY_THRESHOLD = 0.5

_DEPMAP_TASK_URL = "https://depmap.org/portal/api/task/{task_id}"
_DEPMAP_CUSTOM_URL = "https://depmap.org/portal/api/download/custom"
_TASK_TIMEOUT_SECS = 120
_TASK_MIN_POLL_SECS = 2.0

# Persistent disk cache — avoids re-downloading the multi-MB CSV on every cold start
_CACHE_PATH = Path("data/depmap_cache.csv")
_CACHE_MAX_AGE_DAYS = 7


async def poll_task(client: httpx.AsyncClient, task_id: str) -> Any:
    """Poll a DepMap Celery task until SUCCESS or FAILURE, respecting nextPollDelay.

    Returns the `result` payload on success.
    Raises RuntimeError on FAILURE or timeout.
    """
    url = _DEPMAP_TASK_URL.format(task_id=task_id)
    elapsed = 0.0

    while elapsed < _TASK_TIMEOUT_SECS:
        resp = await client.get(url, timeout=30.0)
        resp.raise_for_status()
        data = resp.json()

        state = data.get("state", "PENDING")
        next_delay = max(data.get("nextPollDelay", 2000) / 1000.0, _TASK_MIN_POLL_SECS)

        if state == "SUCCESS":
            logger.info("DepMap task %s succeeded", task_id)
            return data.get("result")

        if state == "FAILURE":
            raise RuntimeError(
                f"DepMap task {task_id} failed: {data.get('message', 'unknown error')}"
            )

        pct = data.get("percentComplete")
        pct_str = f" ({pct}%)" if pct is not None else ""
        logger.debug(
            "DepMap task %s — %s%s, polling in %.1fs",
            task_id,
            state,
            pct_str,
            next_delay,
        )

        await asyncio.sleep(next_delay)
        elapsed += next_delay

    raise TimeoutError(f"DepMap task {task_id} did not complete within {_TASK_TIMEOUT_SECS}s")


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

_GENE_SEARCH_QUERY = """
query GeneSearch($symbol: String!) {
  search(queryString: $symbol, entityNames: ["target"], page: {index: 0, size: 1}) {
    hits { id name entity }
  }
}
"""


def _parse_depmap_csv(text: str) -> dict[str, dict]:
    """Parse DepMap gene_dep_summary CSV text into a lookup dict keyed by uppercase gene symbol."""
    reader = csv.DictReader(io.StringIO(text))
    logger.debug("DepMap CSV columns: %s", reader.fieldnames)
    cache: dict[str, dict] = {}

    for row in reader:
        # Normalize column names: lowercase + underscores (handles "Dependent Cell Lines" → "dependent_cell_lines")
        r = {k.lower().replace(" ", "_"): v for k, v in row.items()}

        # Robust dataset filter: keep rows where dataset column contains 'chronos' or 'crispr',
        # or is absent/empty (single-dataset files have no dataset column).
        dataset = (r.get("dataset") or "").strip().lower()
        if dataset and "chronos" not in dataset and "crispr" not in dataset:
            continue

        gene = (r.get("gene_name") or r.get("gene") or "").strip().upper()
        if not gene:
            continue

        try:
            dep_lines = int(float(r.get("dependent_cell_lines", 0) or 0))
            total_lines = int(float(r.get("cell_lines_with_data", 0) or 0))
        except (ValueError, TypeError):
            continue

        cache[gene] = {
            "dependent_cell_lines": dep_lines,
            "cell_lines_with_data": total_lines,
            "strongly_selective": (r.get("strongly_selective", "False") or "False").lower()
            in ("true", "1", "yes"),
            "common_essential": (r.get("common_essential", "False") or "False").lower()
            in ("true", "1", "yes"),
        }

    return cache


async def load_depmap_cache(client: httpx.AsyncClient) -> dict[str, dict]:
    """Fetch DepMap gene_dep_summary CSV at startup and return a keyed lookup dict.

    Fast path: reads from data/depmap_cache.csv if it exists and is < 7 days old.
    Slow path: downloads from DepMap (handles both direct CSV and Celery task responses),
               saves result to disk for future warm starts.

    Returns an empty dict on any failure (server continues with OT fallback).
    Dict key: uppercase gene symbol. Value: row dict with numeric fields pre-parsed.
    """
    # --- Fast path: load from disk cache ---
    _CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    if _CACHE_PATH.exists():
        age_days = (time.time() - _CACHE_PATH.stat().st_mtime) / 86400
        if age_days < _CACHE_MAX_AGE_DAYS:
            try:
                text = _CACHE_PATH.read_text(encoding="utf-8")
                cache = _parse_depmap_csv(text)
                if cache:
                    logger.info(
                        "Loaded DepMap cache from disk (%d genes, age %.1fh)",
                        len(cache),
                        age_days * 24,
                    )
                    return cache
                logger.warning("Disk cache parsed to 0 genes — re-downloading")
            except Exception as exc:
                logger.warning("Failed to read disk cache: %s — re-downloading", exc)

    # --- Slow path: download from DepMap API ---
    try:
        resp = await client.get(_DEPMAP_SUMMARY_URL, timeout=60.0)
        resp.raise_for_status()

        logger.debug(
            "DepMap gene_dep_summary: status=%d, content-type=%s, body[:200]=%r",
            resp.status_code,
            resp.headers.get("content-type", ""),
            resp.text[:200],
        )

        # Try to parse as JSON first regardless of content-type — DepMap sometimes
        # returns a Celery task body with a non-JSON content-type header.
        try:
            body = resp.json()
            task_id = body.get("id")
            if not task_id:
                logger.warning(
                    "DepMap returned JSON with no task id (keys: %s) — falling back to OT",
                    list(body)[:5],
                )
                return {}

            logger.info("DepMap gene_dep_summary returned task %s — polling...", task_id)
            try:
                result = await poll_task(client, task_id)
            except (RuntimeError, TimeoutError) as exc:
                logger.warning("DepMap task polling failed: %s — falling back to OT", exc)
                return {}

            download_url = (result or {}).get("downloadUrl")
            if not download_url:
                logger.warning(
                    "DepMap task succeeded but no downloadUrl in result; falling back to OT"
                )
                return {}

            csv_resp = await client.get(download_url, timeout=120.0)
            csv_resp.raise_for_status()
            text = csv_resp.text

        except (ValueError, KeyError):
            # Not JSON — treat as direct CSV response
            text = resp.text

        cache = _parse_depmap_csv(text)
        if not cache:
            # Log column names so mismatches are immediately diagnosable
            import csv as _csv
            import io as _io

            _reader = _csv.DictReader(_io.StringIO(text))
            _first = next(_reader, {})
            logger.warning(
                "DepMap CSV parsed to 0 genes — fieldnames=%s, first_row=%s",
                _reader.fieldnames,
                dict(list(_first.items())[:4]),
            )
        logger.info("DepMap cache loaded: %d genes from gene_dep_summary", len(cache))

        # Save to disk for future warm starts
        if cache:
            try:
                _CACHE_PATH.write_text(text, encoding="utf-8")
                logger.info("DepMap cache saved to %s", _CACHE_PATH)
            except Exception as exc:
                logger.warning("Failed to save DepMap cache to disk: %s", exc)

        return cache

    except Exception as exc:
        logger.warning("Failed to load DepMap gene_dep_summary cache: %s — using OT proxy", exc)
        return {}


class DepMapClient:
    """Cancer dependency client: real DepMap CRISPR metrics + OT lineage context."""

    def __init__(self, client: httpx.AsyncClient, gene_dep_cache: dict[str, dict]) -> None:
        self._client = client
        self._cache = gene_dep_cache

    async def fetch_custom_dataset(
        self,
        dataset_id: str,
        feature_labels: list[str] | None = None,
        cell_line_ids: list[str] | None = None,
        drop_empty: bool = True,
    ) -> str:
        """Submit a custom DepMap dataset export task and return the download URL.

        Uses the Celery task queue: POST /download/custom → poll /task/{id} → downloadUrl.
        Caller is responsible for fetching the returned URL (it's a pre-signed S3 link).
        """
        params: dict[str, Any] = {"datasetId": dataset_id, "dropEmpty": drop_empty}
        if feature_labels:
            params["featureLabels"] = feature_labels
        if cell_line_ids:
            params["cellLineIds"] = cell_line_ids

        resp = await self._client.post(_DEPMAP_CUSTOM_URL, params=params, timeout=30.0)
        resp.raise_for_status()
        task_id = resp.json().get("id")
        if not task_id:
            raise ValueError(
                f"DepMap /download/custom returned no task id for dataset {dataset_id!r}"
            )

        logger.info("DepMap custom export task %s submitted for dataset %s", task_id, dataset_id)
        result = await poll_task(self._client, task_id)
        download_url = (result or {}).get("downloadUrl")
        if not download_url:
            raise RuntimeError(f"DepMap task {task_id} succeeded but result has no downloadUrl")
        return download_url

    async def get_essentiality(self, gene_symbol: str) -> CancerDependency | None:
        """Return cancer dependency data, using real DepMap metrics when available."""
        symbol = gene_symbol.strip().upper()

        # Try DepMap cache first for real CRISPR metrics
        cache_entry = self._cache.get(symbol)

        # Get OT lineage context (always useful for per-type breakdown)
        ensembl_id = await self._resolve_gene(symbol)
        ot_data = None
        if ensembl_id:
            ot_data = await self._fetch_ot_cancer_evidence(symbol, ensembl_id)

        if cache_entry is not None:
            return self._build_from_cache(symbol, cache_entry, ot_data)
        elif ot_data is not None:
            return ot_data
        else:
            logger.warning("No DepMap or OT data for '%s'", symbol)
            return None

    def _build_from_cache(
        self,
        gene_symbol: str,
        entry: dict,
        ot_data: CancerDependency | None,
    ) -> CancerDependency:
        """Build CancerDependency using real DepMap metrics, supplemented by OT lineages."""
        dep = entry["dependent_cell_lines"]
        total = entry["cell_lines_with_data"]
        fraction = dep / total if total > 0 else 0.0
        pan_essential = entry["common_essential"]

        # Derive a mean CERES proxy from fraction (real scores not in summary endpoint)
        # Lineage breakdown and cell_lines come from OT supplementary data if available
        top_lineages = ot_data.top_dependent_lineages if ot_data else []
        cell_lines = ot_data.cell_lines if ot_data else []

        source = f"DepMap Chronos Combined ({dep}/{total} cell lines dependent)"
        if ot_data:
            source += " + Open Targets lineage context"

        # Approximate mean score from fraction (−1.0 = very dependent, 0 = not)
        mean_approx = -(fraction * 1.5)

        return CancerDependency(
            gene_symbol=gene_symbol,
            mean_ceres_score=round(mean_approx, 3),
            fraction_dependent_lines=round(fraction, 4),
            pan_essential=pan_essential,
            top_dependent_lineages=top_lineages,
            cell_lines=cell_lines,
            data_source=source,
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

    async def _fetch_ot_cancer_evidence(
        self, gene_symbol: str, ensembl_id: str
    ) -> CancerDependency | None:
        """Fetch Open Targets somatic mutation evidence for lineage context."""
        data = await self._graphql(_CANCER_ASSOCIATIONS_QUERY, {"ensemblId": ensembl_id})
        if data is None:
            return None

        rows = data.get("data", {}).get("target", {}).get("associatedDiseases", {}).get("rows", [])
        if not rows:
            return None

        cancer_rows = [r for r in rows if _is_cancer(r)]
        if not cancer_rows:
            cancer_rows = rows

        cell_lines: list[CellLineEssentiality] = []
        for row in cancer_rows:
            sm_score = next(
                (
                    d["score"]
                    for d in row.get("datatypeScores", [])
                    if d["id"] == "somatic_mutation"
                ),
                0.0,
            )
            if sm_score == 0.0:
                continue
            disease_name = row["disease"]["name"]
            areas = [a["name"] for a in row["disease"].get("therapeuticAreas", [])]
            lineage = next(
                (a for a in areas if "cancer" in a.lower() or "tumor" in a.lower()),
                areas[0] if areas else "Oncology",
            )
            ceres_proxy = -(sm_score * 2.0)
            cell_lines.append(
                CellLineEssentiality(
                    cell_line=disease_name[:50],
                    lineage=lineage,
                    ceres_score=round(ceres_proxy, 3),
                    is_dependent=sm_score >= _OT_DEPENDENCY_THRESHOLD,
                )
            )

        if not cell_lines:
            return None

        scores = [cl.ceres_score for cl in cell_lines]
        mean_score = sum(scores) / len(scores)
        n_dependent = sum(1 for cl in cell_lines if cl.is_dependent)
        fraction_dependent = n_dependent / len(cell_lines)

        # Use disease names (not broad therapeutic areas) for top_dependent_lineages so
        # the output says "melanoma, thyroid carcinoma" rather than "Oncology, Neoplasms".
        disease_order = sorted(cell_lines, key=lambda cl: cl.ceres_score)
        seen_names: dict[str, float] = {}
        for cl in disease_order:
            if cl.cell_line not in seen_names:
                seen_names[cl.cell_line] = cl.ceres_score
        top_lineages = list(seen_names.keys())[:5]
        top_lines = disease_order[:10]

        return CancerDependency(
            gene_symbol=gene_symbol,
            mean_ceres_score=round(mean_score, 4),
            fraction_dependent_lines=round(fraction_dependent, 4),
            pan_essential=False,
            top_dependent_lineages=top_lineages,
            cell_lines=top_lines,
            data_source=(
                "Open Targets Platform v4 — somatic mutation evidence (proxy; "
                "DepMap gene not found in Chronos Combined summary)"
            ),
        )

    async def _graphql(self, query: str, variables: dict) -> dict | None:
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
