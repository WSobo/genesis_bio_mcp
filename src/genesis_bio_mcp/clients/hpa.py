"""Human Protein Atlas (HPA) client.

HPA exposes a per-gene download endpoint that returns a minimal JSON record:

    GET https://www.proteinatlas.org/api/search_download.php?search={symbol}
        &format=json&columns=g,gs,eg,pc,rnatsm,rnats,scml,scl,prognostic_cancer

Columns:
- g       — Ensembl gene ID
- gs      — gene symbol
- eg      — approved synonyms
- pc      — protein class(es); carries the membrane/secreted classification
            ('Predicted membrane proteins', 'Predicted secreted proteins', …)
- rnatsm  — RNA tissue specificity category (e.g. 'Tissue enriched')
- rnats   — RNA tissue specificity score (float)
- scml    — subcellular main location
- scl     — subcellular additional locations
- prognostic_cancer — multiple columns, one per cancer indication

Responses are slow (~2 s/gene) so we disk-cache per symbol with a 7-day TTL,
matching the SAbDab pattern.
"""

from __future__ import annotations

import asyncio
import gzip
import json
import logging
import re
import time
from pathlib import Path

import httpx

from genesis_bio_mcp.config.settings import settings
from genesis_bio_mcp.models import HPAExpression, HPAPathologyData, ProteinAtlasReport

logger = logging.getLogger(__name__)

_HPA_URL = "https://www.proteinatlas.org/api/search_download.php"
_HPA_COLUMNS = "g,gs,eg,pc,rnatsm,rnats,scml,scl,prognostic_cancer"
_SEMAPHORE = asyncio.Semaphore(2)


class HPAClient:
    """Session + disk-cached HPA client.

    Returns ``None`` on unrecoverable errors. Returns an empty
    :class:`ProteinAtlasReport` (expression=None, pathology=[]) when HPA
    has no matching entry for the gene.
    """

    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client
        self._session_cache: dict[str, ProteinAtlasReport] = {}
        self._disk_cache_path: Path = settings.hpa_cache_path
        self._disk_cache: dict[str, dict] = _load_disk_cache(self._disk_cache_path)

    async def get_report(self, gene_symbol: str) -> ProteinAtlasReport | None:
        symbol = gene_symbol.strip().upper()

        if symbol in self._session_cache:
            logger.debug("HPA session cache hit: %s", symbol)
            return self._session_cache[symbol]

        disk_entry = self._disk_cache.get(symbol)
        if (
            disk_entry
            and time.time() - disk_entry.get("fetched_at", 0) < settings.hpa_cache_ttl_secs
        ):
            try:
                report = ProteinAtlasReport(**disk_entry["report"])
                self._session_cache[symbol] = report
                return report
            except Exception as exc:
                logger.debug("HPA disk cache entry for %s stale or malformed: %s", symbol, exc)

        async with _SEMAPHORE:
            raw = await self._fetch(symbol)
        if raw is None:
            return None

        report = _parse_hpa(raw, symbol)
        self._session_cache[symbol] = report
        if report.expression is not None or report.pathology:
            self._disk_cache[symbol] = {
                "fetched_at": time.time(),
                "report": report.model_dump(),
            }
            _save_disk_cache(self._disk_cache_path, self._disk_cache)
        return report

    async def _fetch(self, symbol: str) -> dict | None:
        params = {
            "search": symbol,
            "format": "json",
            "columns": _HPA_COLUMNS,
        }
        try:
            resp = await self._client.get(_HPA_URL, params=params, timeout=25.0)
            if resp.status_code == 404:
                return None
            resp.raise_for_status()
            # HPA's search_download.php returns a gzip-compressed body
            # (Content-Type: application/gzip) — a file payload, NOT a
            # Content-Encoding transport that httpx would auto-decode. So the
            # raw body starts with the gzip magic bytes (0x1f 0x8b) and
            # resp.json() would fail. Detect and inflate before parsing.
            content = resp.content
            if content[:2] == b"\x1f\x8b":
                content = gzip.decompress(content)
            data = json.loads(content)
        except Exception as exc:
            logger.warning("HPA fetch failed for %s: %s", symbol, exc)
            return None

        # HPA returns a list; the first element matching the exact symbol is
        # the one we want. Free-text search can return near-neighbours.
        if not isinstance(data, list) or not data:
            return {}
        for row in data:
            if str(row.get("Gene") or row.get("gs") or "").upper() == symbol:
                return row
        return data[0]


_PROGNOSTIC_KEY = re.compile(r"^Pathology prognostics - (.+)$", re.IGNORECASE)


def _as_str_list(value: object) -> list[str]:
    """Normalize an HPA multi-value field (list or comma-separated string) to a list.

    HPA's JSON download returns some columns as arrays and others as a single
    comma-separated string; this collapses both into a clean, de-duplicated list.
    """
    if value in (None, ""):
        return []
    if isinstance(value, list):
        items = [str(v).strip() for v in value]
    else:
        items = [s.strip() for s in str(value).split(",")]
    seen: set[str] = set()
    return [s for s in items if s and not (s in seen or seen.add(s))]


def _parse_hpa(row: dict, symbol: str) -> ProteinAtlasReport:
    """Convert one HPA row into a :class:`ProteinAtlasReport`."""
    if not row:
        return ProteinAtlasReport(gene_symbol=symbol)

    ensembl_id = row.get("Ensembl") or row.get("g")
    specificity_cat = row.get("RNA tissue specificity") or row.get("rnatsm")
    specificity_score = row.get("RNA tissue specificity score") or row.get("rnats")
    try:
        spec_score_f = float(specificity_score) if specificity_score not in (None, "") else None
    except (TypeError, ValueError):
        spec_score_f = None

    # Subcellular locations: HPA's JSON download returns these as arrays (older
    # shapes used comma-separated strings). _as_str_list normalizes both; merge
    # main + additional and de-duplicate while preserving order. (str() on a list
    # would otherwise leak a Python list repr into the output.)
    subcellular_raw = _as_str_list(row.get("Subcellular main location") or row.get("scml"))
    subcellular_raw += _as_str_list(row.get("Subcellular location") or row.get("scl"))
    seen: set[str] = set()
    subcellular = [s for s in subcellular_raw if not (s in seen or seen.add(s))]

    # Enhanced tissues: HPA returns "RNA tissue specific nTPM" as a JSON object
    # ({tissue: nTPM}); render "tissue (nTPM)". Fall back to list/string shapes.
    enhanced_raw = row.get("RNA tissue specific nTPM") or row.get("Tissue expression cluster")
    if isinstance(enhanced_raw, dict):
        enhanced_tissues = [f"{k} ({v})" for k, v in enhanced_raw.items()]
    elif isinstance(enhanced_raw, list):
        enhanced_tissues = [str(t).strip() for t in enhanced_raw if str(t).strip()]
    else:
        enhanced_tissues = [t.strip() for t in str(enhanced_raw or "").split(";") if t.strip()]

    # Protein class — HPA returns this as a list (JSON) or a comma-separated
    # string. It carries the membrane/secreted classification, the strongest
    # single signal for whether a target is reachable by antibodies/biologics.
    protein_classes = _as_str_list(row.get("Protein class") or row.get("pc"))
    pc_lower = " ".join(protein_classes).lower()
    is_membrane = "membrane" in pc_lower
    is_secreted = "secreted" in pc_lower

    expression = HPAExpression(
        gene_symbol=symbol,
        ensembl_id=str(ensembl_id) if ensembl_id else None,
        rna_tissue_specificity_category=str(specificity_cat) if specificity_cat else None,
        rna_tissue_specificity_score=spec_score_f,
        enhanced_tissues=enhanced_tissues,
        subcellular_locations=subcellular,
        protein_classes=protein_classes,
        is_membrane=is_membrane,
        is_secreted=is_secreted,
    )

    pathology: list[HPAPathologyData] = []
    for key, value in row.items():
        m = _PROGNOSTIC_KEY.match(str(key))
        if not m or not value:
            continue
        text = str(value).strip()
        if not text or text.lower() in ("none", "not significant"):
            continue
        prognostic = None
        low = text.lower()
        # Check "unfavorable" first — "favorable" is a substring of it.
        if "unfavorable" in low or "unfavourable" in low:
            prognostic = "Unfavorable"
        elif "favorable" in low or "favourable" in low:
            prognostic = "Favorable"
        pathology.append(
            HPAPathologyData(
                cancer_type=m.group(1).strip(),
                prognostic_outcome=prognostic,
                staining_intensity=None,
            )
        )

    return ProteinAtlasReport(gene_symbol=symbol, expression=expression, pathology=pathology)


def _load_disk_cache(path: Path) -> dict[str, dict]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except Exception as exc:
        logger.warning("HPA disk cache unreadable at %s: %s", path, exc)
        return {}


def _save_disk_cache(path: Path, cache: dict[str, dict]) -> None:
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(cache))
    except Exception as exc:
        logger.warning("HPA disk cache write failed at %s: %s", path, exc)
