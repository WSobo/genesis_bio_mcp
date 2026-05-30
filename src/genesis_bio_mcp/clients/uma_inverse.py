"""UMA-Inverse — deployed inverse-folding service client (structure → sequence).

Wraps the REST service (default https://wsobo-uma-inverse.hf.space, configurable via
``GENESIS_UMA_API_URL``):

  - ``POST /design``  — redesign a backbone's sequence.
  - ``POST /score``   — score a sequence against a structure (→ candidate mutations).

This is the stack's first *served-model* client: it consumes a deployed, monitored ML
endpoint rather than a database. The structure can be supplied as raw PDB text or as a URL
(``resolve_pdb`` fetches it) — so an agent can pass the AlphaFold model URL that
``get_protein_structure`` returns. Never raises to the caller; returns ``None`` on any
failure (unreachable service, over the residue cap → 413, bad body → 422, timeout → 504).
"""

from __future__ import annotations

import asyncio
import logging

import httpx

from genesis_bio_mcp.config.settings import settings
from genesis_bio_mcp.models import ScoredPosition, SequenceDesign, StructureScore

logger = logging.getLogger(__name__)

# A single CPU Space serves one request at a time; keep concurrency low.
_SEMAPHORE = asyncio.Semaphore(2)


class UMAInverseClient:
    """Client for the deployed UMA-Inverse inverse-folding REST service."""

    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client

    @property
    def _base(self) -> str:
        return settings.uma_api_url.rstrip("/")

    async def resolve_pdb(self, structure: str) -> str | None:
        """Return PDB text. If *structure* is a URL, fetch it; otherwise return it as-is."""
        s = (structure or "").strip()
        if not s:
            return None
        if s.startswith(("http://", "https://")):
            try:
                resp = await self._client.get(
                    s, timeout=settings.uma_timeout_secs, follow_redirects=True
                )
                resp.raise_for_status()
                return resp.text
            except Exception as exc:
                logger.warning("UMA: failed to fetch PDB from %s: %s", s, exc)
                return None
        return s

    async def design(
        self,
        pdb: str,
        *,
        ligand: str | None = None,
        temperature: float = 0.1,
        n_samples: int = 1,
    ) -> SequenceDesign | None:
        """Redesign a backbone's sequence via ``POST /design``."""
        data = await self._post(
            "/design",
            {"pdb": pdb, "ligand": ligand, "temperature": temperature, "n_samples": n_samples},
        )
        if data is None:
            return None
        return SequenceDesign(
            n_residues=int(data.get("n_residues") or 0),
            mean_confidence=data.get("mean_confidence"),
            sequences=data.get("sequences") or [],
            per_residue_confidence=data.get("per_residue_confidence") or [],
            inference_ms=data.get("inference_ms"),
            request_id=data.get("request_id"),
        )

    async def score(
        self,
        pdb: str,
        *,
        sequence: str | None = None,
        mode: str = "autoregressive",
    ) -> StructureScore | None:
        """Score a sequence against a structure via ``POST /score``."""
        data = await self._post("/score", {"pdb": pdb, "sequence": sequence, "mode": mode})
        if data is None:
            return None
        positions = [
            ScoredPosition(
                position=int(p.get("position") or 0),
                residue_id=p.get("residue_id"),
                aa=p.get("aa"),
                log_prob=p.get("log_prob"),
                prob=p.get("prob"),
                top_aa=p.get("top_aa"),
                top_prob=p.get("top_prob"),
            )
            for p in (data.get("positions") or [])
        ]
        return StructureScore(
            n_residues=int(data.get("n_residues") or 0),
            perplexity=data.get("perplexity"),
            recovery=data.get("recovery"),
            mean_log_prob=data.get("mean_log_prob"),
            mode=data.get("mode"),
            sequence_scored=data.get("sequence_scored"),
            positions=positions,
            inference_ms=data.get("inference_ms"),
            request_id=data.get("request_id"),
        )

    async def _post(self, path: str, payload: dict) -> dict | None:
        url = f"{self._base}{path}"
        try:
            async with _SEMAPHORE:
                resp = await self._client.post(url, json=payload, timeout=settings.uma_timeout_secs)
            if resp.status_code != 200:
                logger.warning("UMA %s returned %s: %s", path, resp.status_code, resp.text[:200])
                return None
            return resp.json()
        except Exception as exc:
            logger.warning("UMA %s request failed: %s", path, exc)
            return None
