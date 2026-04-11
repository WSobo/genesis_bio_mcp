"""EFO ontology resolver via EBI OLS4 API.

Converts a free-text trait query (including colloquial terms) to EFO ontology
terms, enabling ontology-backed matching against GWAS Catalog associations.

Why this beats a synonym dict:
- EFO is maintained by EBI ontologists who track GWAS Catalog vocabulary.
- Synonyms in EFO cover clinical terms, colloquialisms, and non-English variants.
- New indications work automatically without touching code.
- Hierarchy: matching a parent EFO term captures child terms in GWAS results.

Cache strategy: 7-day disk cache (EFO releases ~quarterly; daily invalidation
is wasteful). Session cache avoids repeated OLS calls within one MCP session.
"""

from __future__ import annotations

import json
import logging
import time
import unicodedata
from dataclasses import dataclass, field
from pathlib import Path

import httpx

logger = logging.getLogger(__name__)

_OLS_URL = "https://www.ebi.ac.uk/ols4/api/search"
_EFO_CACHE_PATH = Path("data/efo_cache.json")
_EFO_CACHE_TTL_SECS = 604800  # 7 days


def _normalize(text: str) -> str:
    return unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii").lower()


@dataclass
class EFOTerm:
    uri: str  # "http://www.ebi.ac.uk/efo/EFO_0001073"
    label: str  # "obesity"
    synonyms: list[str] = field(default_factory=list)


def _load_efo_cache() -> dict[str, dict]:
    try:
        if _EFO_CACHE_PATH.exists():
            return json.loads(_EFO_CACHE_PATH.read_text())
    except Exception as exc:
        logger.warning("Failed to load EFO cache: %s", repr(exc))
    return {}


def _parse_docs(docs: list[dict]) -> list[EFOTerm]:
    terms: list[EFOTerm] = []
    for doc in docs:
        uri = doc.get("iri", "")
        label = doc.get("label", "")
        if not uri or not label:
            continue
        # OLS4 returns synonym as list[str] or absent
        raw_syns = doc.get("synonym", []) or []
        synonyms = [s for s in raw_syns if isinstance(s, str) and s]
        terms.append(EFOTerm(uri=uri, label=label, synonyms=synonyms))
    return terms


class EFOResolver:
    """Resolve free-text trait queries to EFO terms via OLS4.

    Call path: session cache → disk cache → OLS4 API.
    Never raises — returns [] on any failure so callers fall back to synonyms.

    Args:
        client: Shared httpx.AsyncClient from the server lifespan.
        cache_path: Override the disk cache location. Pass ``None`` to disable
            disk caching entirely (useful for tests or ephemeral deployments).
    """

    def __init__(
        self,
        client: httpx.AsyncClient,
        *,
        cache_path: Path | None = _EFO_CACHE_PATH,
    ) -> None:
        self._client = client
        self._cache_path = cache_path
        self._session_cache: dict[str, list[EFOTerm]] = {}
        self._disk_cache: dict[str, dict] = _load_efo_cache() if cache_path else {}

    async def resolve(self, trait: str) -> list[EFOTerm]:
        """Return up to 5 best-matching EFO terms for the trait query."""
        key = _normalize(trait)

        # 1. Session cache
        if key in self._session_cache:
            return self._session_cache[key]

        # 2. Disk cache (7-day TTL)
        entry = self._disk_cache.get(key)
        if entry and time.time() - entry.get("fetched_at", 0) < _EFO_CACHE_TTL_SECS:
            try:
                terms = [EFOTerm(**t) for t in entry["terms"]]
                self._session_cache[key] = terms
                return terms
            except Exception:
                pass

        # 3. OLS4 API
        terms = await self._fetch_from_ols(trait, key)
        self._session_cache[key] = terms
        self._write_disk_cache(key, terms)
        return terms

    async def _fetch_from_ols(self, trait: str, key: str) -> list[EFOTerm]:
        try:
            resp = await self._client.get(
                _OLS_URL,
                params={
                    "q": trait,
                    "ontology": "efo",
                    "type": "class",
                    "rows": 5,
                    "fieldList": "iri,label,synonym",
                },
                timeout=8.0,
            )
            resp.raise_for_status()
            docs = resp.json().get("response", {}).get("docs", [])
            terms = _parse_docs(docs)
            logger.info(
                "OLS4 resolved '%s' → %d EFO terms: %s",
                trait,
                len(terms),
                [t.label for t in terms],
            )
            return terms
        except Exception as exc:
            logger.warning("OLS4 lookup failed for '%s': %s", trait, repr(exc))
            return []

    def _write_disk_cache(self, key: str, terms: list[EFOTerm]) -> None:
        if not self._cache_path:
            return
        self._disk_cache[key] = {
            "terms": [{"uri": t.uri, "label": t.label, "synonyms": t.synonyms} for t in terms],
            "fetched_at": time.time(),
        }
        try:
            self._cache_path.parent.mkdir(parents=True, exist_ok=True)
            self._cache_path.write_text(json.dumps(self._disk_cache, indent=2, default=str))
        except Exception as exc:
            logger.warning("Failed to write EFO cache: %s", repr(exc))
