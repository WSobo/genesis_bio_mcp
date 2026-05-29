"""Foldseek structural-similarity search-server client.

Wraps the public Foldseek web search server (https://search.foldseek.com).
It is an asynchronous job queue:

1. ``POST /api/ticket``        — submit a structure (multipart ``q`` file +
   ``database[]`` + ``mode``); returns ``{"id", "status"}``.
2. ``GET  /api/ticket/{id}``   — poll status: PENDING | RUNNING | COMPLETE | ERROR.
3. ``GET  /api/result/{id}/0`` — fetch the alignment JSON once COMPLETE.

The query structure is the gene's AlphaFold model PDB, fetched via the
injected :class:`AlphaFoldClient`. Foldseek is a shared academic resource, so
the module-level semaphore is pinned to 1 concurrent search.
"""

from __future__ import annotations

import asyncio
import logging
import re

import httpx

from genesis_bio_mcp.clients.alphafold import AlphaFoldClient
from genesis_bio_mcp.config.settings import settings
from genesis_bio_mcp.models import StructuralHomolog, StructuralHomologs

logger = logging.getLogger(__name__)

_BASE_URL = "https://search.foldseek.com/api"
# Foldseek's default 3Di+AA mode — fast and appropriate for fold-level homology.
_SEARCH_MODE = "3diaa"
# Public shared server — never run two searches at once.
_SEMAPHORE = asyncio.Semaphore(1)

# AlphaFold-DB target ids look like ``AF-P15056-F1-model_v4``; pull the accession.
_AFDB_ACCESSION_RE = re.compile(r"^AF-([A-Z0-9]+)-F\d+")


class FoldseekClient:
    """Structural-homolog search via the Foldseek web server."""

    def __init__(self, client: httpx.AsyncClient, *, alphafold: AlphaFoldClient) -> None:
        self._client = client
        self._alphafold = alphafold
        # Session cache keyed by (accession, database). Only successes are stored.
        self._cache: dict[tuple[str, str], StructuralHomologs] = {}

    async def search(
        self,
        gene_symbol: str,
        uniprot_accession: str | None,
        *,
        database: str | None = None,
        max_hits: int | None = None,
    ) -> StructuralHomologs | None:
        """Return structurally-similar proteins for a gene's AlphaFold model.

        Resolves the gene's AlphaFold model PDB, submits it to Foldseek, polls
        the ticket to completion, and parses the top hits. Returns ``None`` on
        any failure (no model, submission/poll/parse error, timeout) — never
        raises to the caller.
        """
        if not uniprot_accession:
            logger.warning("Foldseek: no UniProt accession for %s", gene_symbol)
            return None

        db = database or settings.foldseek_database
        limit = max_hits if max_hits is not None else settings.foldseek_max_hits

        cache_key = (uniprot_accession, db)
        if cache_key in self._cache:
            logger.debug("Foldseek cache hit: %s/%s", uniprot_accession, db)
            cached = self._cache[cache_key]
            # Re-slice to the requested limit without re-querying.
            return cached.model_copy(update={"hits": cached.hits[:limit]})

        pdb_text = await self._alphafold.get_model_pdb(uniprot_accession)
        if not pdb_text:
            logger.info("Foldseek: no AlphaFold model for %s (%s)", gene_symbol, uniprot_accession)
            return None

        async with _SEMAPHORE:
            ticket = await self._submit(pdb_text, db)
            if ticket is None:
                return None
            if not await self._poll(ticket):
                return None
            body = await self._fetch_results(ticket)

        if body is None:
            return None

        hits = _parse_results(body, uniprot_accession, limit)
        result = StructuralHomologs(
            gene_symbol=gene_symbol,
            uniprot_accession=uniprot_accession,
            database=db,
            hits=hits,
        )
        # Cache the full (unsliced-by-this-call) result; store under the db key.
        self._cache[cache_key] = result
        return result

    async def _submit(self, pdb_text: str, database: str) -> str | None:
        """Submit a structure search ticket; return the ticket id."""
        try:
            resp = await self._client.post(
                f"{_BASE_URL}/ticket",
                files={"q": ("query.pdb", pdb_text.encode(), "application/octet-stream")},
                data={"database[]": database, "mode": _SEARCH_MODE},
                timeout=30.0,
            )
            resp.raise_for_status()
            ticket = resp.json().get("id")
            if not ticket:
                logger.warning("Foldseek submit returned no ticket id")
                return None
            return str(ticket)
        except Exception as exc:
            logger.warning("Foldseek submit failed: %s", exc)
            return None

    async def _poll(self, ticket: str) -> bool:
        """Poll a ticket until COMPLETE. Return True on success, False otherwise."""
        elapsed = 0.0
        interval = settings.foldseek_poll_interval_secs
        while elapsed < settings.foldseek_poll_timeout_secs:
            try:
                resp = await self._client.get(f"{_BASE_URL}/ticket/{ticket}", timeout=20.0)
                resp.raise_for_status()
                status = resp.json().get("status")
            except Exception as exc:
                logger.warning("Foldseek poll failed for %s: %s", ticket, exc)
                return False
            if status == "COMPLETE":
                return True
            if status == "ERROR":
                logger.warning("Foldseek job %s reported ERROR", ticket)
                return False
            await asyncio.sleep(interval)
            elapsed += interval
        logger.warning(
            "Foldseek job %s timed out after %ss", ticket, settings.foldseek_poll_timeout_secs
        )
        return False

    async def _fetch_results(self, ticket: str) -> dict | None:
        """Fetch the result JSON for a completed ticket."""
        try:
            resp = await self._client.get(f"{_BASE_URL}/result/{ticket}/0", timeout=30.0)
            resp.raise_for_status()
            return resp.json()
        except Exception as exc:
            logger.warning("Foldseek result fetch failed for %s: %s", ticket, exc)
            return None


def _split_target(target: str) -> tuple[str, str | None]:
    """Split a Foldseek target header ``"<id> <description>"`` into (id, description)."""
    target = (target or "").strip()
    if not target:
        return "", None
    parts = target.split(None, 1)
    target_id = parts[0]
    description = parts[1].strip() if len(parts) > 1 else None
    return target_id, description


def _safe_float(value: object) -> float | None:
    try:
        return float(value) if value is not None else None  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def _parse_results(body: dict, query_accession: str, max_hits: int) -> list[StructuralHomolog]:
    """Parse the Foldseek result JSON into a sorted, self-excluded hit list.

    ``results`` is a list of per-database blocks; each block's ``alignments`` is
    a list-of-lists (one inner list per submitted query — we submit a single
    query, so we read the first). The query's own AlphaFold model is excluded.
    """
    hits: list[StructuralHomolog] = []
    for block in body.get("results") or []:
        db = block.get("db", "")
        alignments = block.get("alignments") or []
        # Flatten the per-query nesting (alignments[0] is this query's hit list).
        flat = alignments[0] if alignments and isinstance(alignments[0], list) else alignments
        for a in flat or []:
            target = a.get("target", "")
            # Drop the self-hit: the query's own AlphaFold model.
            if query_accession and f"AF-{query_accession}-" in target:
                continue
            target_id, description = _split_target(target)
            accession_match = _AFDB_ACCESSION_RE.match(target_id)
            hits.append(
                StructuralHomolog(
                    target_id=target_id,
                    description=description,
                    database=db,
                    uniprot_accession=accession_match.group(1) if accession_match else None,
                    evalue=_safe_float(a.get("eval")),
                    bit_score=_safe_float(a.get("score")),
                    probability=_safe_float(a.get("prob")),
                    seq_identity=_safe_float(a.get("seqId")),
                    aln_length=a.get("alnLength"),
                )
            )
    hits.sort(key=lambda h: h.evalue if h.evalue is not None else float("inf"))
    return hits[:max_hits]
