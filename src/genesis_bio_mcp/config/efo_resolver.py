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

import asyncio
import json
import logging
import time
import unicodedata
from dataclasses import dataclass, field
from pathlib import Path

import httpx

from genesis_bio_mcp.config.settings import settings

logger = logging.getLogger(__name__)

_OLS_URL = "https://www.ebi.ac.uk/ols4/api/search"
# Hierarchy expansion bounds. Descendants are large-ish (disease subtype clades) but
# capped to avoid pathological ontology nodes. Ancestors are capped tightly so that
# "polycythemia vera" doesn't expand all the way up to "disease" and catch every
# disease-tagged study.
_MAX_DESCENDANTS = 500
_MAX_ANCESTORS = 10

# Bump when the cached entry schema changes so stale entries are invalidated
# on first read rather than served without required fields.
# v2 (0.2.5): adds ``related_uris`` populated by hierarchy expansion. Without
# the bump, entries written by v0.2.3 load with an empty list (dataclass
# default) and the expansion is never re-triggered.
_CACHE_SCHEMA_VERSION = 2
# Default cache path read from settings at class definition time so it can be
# used as the default argument to EFOResolver.__init__ below.
_EFO_CACHE_PATH: Path = settings.efo_cache_path


def _normalize(text: str) -> str:
    return unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii").lower()


@dataclass
class EFOTerm:
    uri: str  # "http://www.ebi.ac.uk/efo/EFO_0001073"
    label: str  # "obesity"
    synonyms: list[str] = field(default_factory=list)
    # Hierarchy-expanded URIs: descendants (more specific subtypes) + direct ancestors
    # (one level broader). Populated by EFOResolver.resolve() via OLS4 structural
    # filters. Feeds the URI-set match in filter_by_trait so queries like
    # "polycythemia vera" catch JAK2 studies tagged "myeloproliferative neoplasm"
    # (direct parent) and queries like "myeloproliferative disorder" catch studies
    # tagged with specific subtypes.
    related_uris: list[str] = field(default_factory=list)


def _load_efo_cache(cache_path: Path) -> dict[str, dict]:
    try:
        if cache_path.exists():
            return json.loads(cache_path.read_text())
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
        self._disk_cache: dict[str, dict] = _load_efo_cache(cache_path) if cache_path else {}

    async def resolve(self, trait: str) -> list[EFOTerm]:
        """Return up to 5 best-matching EFO terms for the trait query.

        Each returned term is hierarchy-expanded: its ``related_uris`` contains
        URIs of descendants (for broader queries like "myeloproliferative") and
        direct ancestors (for specific queries like "polycythemia vera"), so
        ontology-backed URI matching catches siblings/parents/children of the
        literal query term.
        """
        key = _normalize(trait)

        # 1. Session cache
        if key in self._session_cache:
            return self._session_cache[key]

        # 2. Disk cache (7-day TTL). Reject entries written by earlier schema
        # versions so the hierarchy expansion runs on the next resolve after
        # an upgrade. Without this check, v0.2.3 cache entries load fine
        # (dataclass default fills related_uris with []) and the expansion
        # is never re-triggered, which was the root cause of the v0.2.4
        # "polycythemia vera still misses" regression.
        entry = self._disk_cache.get(key)
        if (
            entry
            and entry.get("schema_version", 0) >= _CACHE_SCHEMA_VERSION
            and time.time() - entry.get("fetched_at", 0) < settings.efo_cache_ttl_secs
        ):
            try:
                terms = [EFOTerm(**t) for t in entry["terms"]]
                self._session_cache[key] = terms
                return terms
            except Exception:
                pass

        # 3. OLS4 API — fetch matching terms, then expand each term's hierarchy
        terms = await self._fetch_from_ols(trait, key)
        if terms:
            await self._expand_hierarchy(terms)
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

    async def _expand_hierarchy(self, terms: list[EFOTerm]) -> None:
        """Populate ``related_uris`` on each term with descendants + direct ancestors.

        Uses OLS4 search with ``allChildrenOf`` / ``ancestorsOf`` filters to fetch
        structural neighbors of each term. Runs all fetches concurrently. Any
        failure silently leaves ``related_uris`` empty — callers already handle
        sparse ontology matches via the synonym fallback in filter_by_trait.
        """
        # Fan out one (desc, anc) fetch pair per term; gather all results in parallel
        coros = []
        for t in terms:
            coros.append(self._fetch_related(t.uri, "allChildrenOf", _MAX_DESCENDANTS))
            coros.append(self._fetch_related(t.uri, "ancestorsOf", _MAX_ANCESTORS))
        results = await asyncio.gather(*coros)

        # Stitch results back onto their owning terms (descendants, ancestors per term)
        for i, t in enumerate(terms):
            desc = results[2 * i]
            anc = results[2 * i + 1]
            merged = set(desc) | set(anc)
            merged.discard(t.uri)
            t.related_uris = sorted(merged)

    async def _fetch_related(self, iri: str, filter_param: str, rows: int) -> list[str]:
        """OLS4 search restricted by ``filter_param`` (``allChildrenOf``/``ancestorsOf``)."""
        try:
            resp = await self._client.get(
                _OLS_URL,
                params={
                    "q": "*",
                    filter_param: iri,
                    "ontology": "efo",
                    "type": "class",
                    "rows": rows,
                    "fieldList": "iri",
                },
                timeout=8.0,
            )
            resp.raise_for_status()
            docs = resp.json().get("response", {}).get("docs", [])
            return [d.get("iri", "") for d in docs if d.get("iri")]
        except Exception as exc:
            logger.warning(
                "OLS4 hierarchy lookup failed (%s for %s): %s", filter_param, iri, repr(exc)
            )
            return []

    def _write_disk_cache(self, key: str, terms: list[EFOTerm]) -> None:
        if not self._cache_path:
            return
        self._disk_cache[key] = {
            "schema_version": _CACHE_SCHEMA_VERSION,
            "terms": [
                {
                    "uri": t.uri,
                    "label": t.label,
                    "synonyms": t.synonyms,
                    "related_uris": t.related_uris,
                }
                for t in terms
            ],
            "fetched_at": time.time(),
        }
        try:
            self._cache_path.parent.mkdir(parents=True, exist_ok=True)
            self._cache_path.write_text(json.dumps(self._disk_cache, indent=2, default=str))
        except Exception as exc:
            logger.warning("Failed to write EFO cache: %s", repr(exc))
