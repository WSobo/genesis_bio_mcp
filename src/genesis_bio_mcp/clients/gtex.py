"""GTEx API client — bulk-RNA median gene expression per tissue.

Wraps https://gtexportal.org/api/v2. No API key required.

GTEx requires the **versioned** GENCODE ID on the expression endpoint
(``ENSG00000169174.11``); the unversioned form silently returns an empty
``data: []`` payload.

The expression dataset (``gtex_v10``) indexes against **GENCODE v39**, but
``/reference/gene`` defaults to GENCODE v26 (the historical default) and
returns the v26-style ID — which the v39-keyed expression endpoint then
fails to find. The reference endpoint accepts ``gencodeVersion`` ∈
{``v19``, ``v26``, ``v39``}; we pin it to ``v39`` so the returned ID
matches what the expression endpoint serves.

Endpoints:
- ``GET /api/v2/reference/gene?geneId={symbol}&gencodeVersion=v39``  — symbol → v39 GENCODE ID
- ``GET /api/v2/expression/medianGeneExpression?gencodeId={id}``     — median TPM per tissue
"""

from __future__ import annotations

import asyncio
import json
import logging
import time
from pathlib import Path

import httpx

from genesis_bio_mcp.clients.ensembl import EnsemblClient
from genesis_bio_mcp.config.settings import settings
from genesis_bio_mcp.models import GTExExpression, TissueExpressionProfile

logger = logging.getLogger(__name__)

_GTEX_BASE = "https://gtexportal.org/api/v2"
_SEMAPHORE = asyncio.Semaphore(3)


class GTExClient:
    """Session + disk-cached GTEx client.

    Disk cache stores the full parsed :class:`TissueExpressionProfile` per gene
    symbol. Keyed by symbol (not GENCODE ID) so repeated lookups with different
    GENCODE versions hit the same cache entry.
    """

    def __init__(self, client: httpx.AsyncClient, *, ensembl: EnsemblClient) -> None:
        self._client = client
        self._ensembl = ensembl
        self._session_cache: dict[str, TissueExpressionProfile] = {}
        self._disk_cache_path: Path = settings.gtex_cache_path
        self._disk_cache: dict[str, dict] = _load_disk_cache(self._disk_cache_path)

    async def get_expression(self, gene_symbol: str) -> TissueExpressionProfile | None:
        """Return median TPM per tissue for *gene_symbol*.

        Returns ``None`` only on unrecoverable errors. Returns a profile with
        an empty ``samples`` list when the gene has no GTEx mapping (e.g.
        non-coding genes, genes retired from GENCODE v10+).
        """
        symbol = gene_symbol.strip().upper()

        if symbol in self._session_cache:
            logger.debug("GTEx session cache hit: %s", symbol)
            return self._session_cache[symbol]

        disk_entry = self._disk_cache.get(symbol)
        if (
            disk_entry
            and time.time() - disk_entry.get("fetched_at", 0) < settings.gtex_cache_ttl_secs
        ):
            try:
                profile = TissueExpressionProfile(**disk_entry["profile"])
                self._session_cache[symbol] = profile
                return profile
            except Exception as exc:
                logger.debug("GTEx disk cache entry for %s stale or malformed: %s", symbol, exc)

        # Resolve symbol → versioned GENCODE ID via GTEx's own reference
        # endpoint. The version suffix (e.g. ``.11``) is mandatory on the
        # expression endpoint; without it GTEx returns an empty data array
        # rather than 404. Ensembl is the secondary fallback for symbols
        # GTEx doesn't index (rare).
        gencode_id = await self._resolve_gencode_id(symbol)
        if gencode_id is None:
            ensembl_gene = await self._ensembl.lookup_gene(symbol)
            if ensembl_gene is not None:
                gencode_id = ensembl_gene.ensembl_id
        if gencode_id is None:
            logger.debug("GTEx: no GENCODE ID for %s", symbol)
            empty = TissueExpressionProfile(gene_symbol=symbol, gencode_id=None)
            self._session_cache[symbol] = empty
            return empty

        async with _SEMAPHORE:
            samples = await self._fetch_expression(gencode_id)

        profile = TissueExpressionProfile(
            gene_symbol=symbol, gencode_id=gencode_id, samples=samples
        )
        self._session_cache[symbol] = profile
        # Persist only successful fetches (samples non-empty) so empty rows
        # are retried on future calls.
        if samples:
            self._disk_cache[symbol] = {
                "fetched_at": time.time(),
                "profile": profile.model_dump(),
            }
            _save_disk_cache(self._disk_cache_path, self._disk_cache)
        return profile

    async def _resolve_gencode_id(self, gene_symbol: str) -> str | None:
        """Map an HGNC symbol to GTEx's pinned versioned GENCODE ID.

        Returns ``None`` if GTEx doesn't index the symbol (common for
        non-coding / pseudogene aliases) or on any network error.
        """
        url = f"{_GTEX_BASE}/reference/gene"
        try:
            async with _SEMAPHORE:
                resp = await self._client.get(
                    url,
                    params={"geneId": gene_symbol, "gencodeVersion": "v39"},
                    timeout=20.0,
                )
            if resp.status_code == 404:
                return None
            resp.raise_for_status()
            data = resp.json()
        except Exception as exc:
            logger.warning("GTEx /reference/gene failed for %s: %s", gene_symbol, exc)
            return None
        rows = data.get("data") or []
        for row in rows:
            # GTEx returns one row per matched symbol; prefer the entry whose
            # primary geneSymbol matches the query exactly to avoid
            # near-miss matches (e.g. "ALB" vs "ALB2") biting us the way the
            # UniProt resolver did.
            if (row.get("geneSymbol") or "").upper() == gene_symbol.strip().upper():
                gencode = row.get("gencodeId")
                if gencode:
                    return str(gencode)
        # Fallback to first row if no exact symbol match — GTEx's matching
        # is symbol-prefix-based and usually returns the canonical entry first.
        if rows:
            return rows[0].get("gencodeId")
        return None

    async def _fetch_expression(self, gencode_id: str) -> list[GTExExpression]:
        """Call /api/v2/expression/medianGeneExpression and parse the rows."""
        url = f"{_GTEX_BASE}/expression/medianGeneExpression"
        params = {"gencodeId": gencode_id}
        try:
            resp = await self._client.get(url, params=params, timeout=25.0)
            if resp.status_code == 404:
                return []
            resp.raise_for_status()
            data = resp.json()
        except Exception as exc:
            logger.warning("GTEx fetch failed for %s: %s", gencode_id, exc)
            return []

        rows = data.get("data") or data.get("medianGeneExpression") or []
        samples: list[GTExExpression] = []
        for row in rows:
            tissue = row.get("tissueSiteDetailId") or row.get("tissue") or ""
            # GTEx returns underscored labels like "Brain_Cortex"; prettify to
            # match HPA's hyphenated form, which tooling expects.
            if tissue:
                tissue = tissue.replace("_", " - ")
            median = row.get("median") or row.get("medianTpm")
            if tissue and median is not None:
                try:
                    samples.append(
                        GTExExpression(
                            tissue=tissue,
                            median_tpm=float(median),
                            sample_count=row.get("sampleCount") or row.get("numSamples"),
                        )
                    )
                except (TypeError, ValueError):
                    continue
        return samples


def _load_disk_cache(path: Path) -> dict[str, dict]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except Exception as exc:
        logger.warning("GTEx disk cache unreadable at %s: %s", path, exc)
        return {}


def _save_disk_cache(path: Path, cache: dict[str, dict]) -> None:
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(cache))
    except Exception as exc:
        logger.warning("GTEx disk cache write failed at %s: %s", path, exc)
