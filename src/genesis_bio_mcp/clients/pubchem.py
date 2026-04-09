"""PubChem REST + NCBI E-utils client for compound bioactivity data.

PubChem's concisebioactivities endpoint is invalid for gene/assay targets,
and most BRAF assays (1456 AIDs) use continuous binding formats with no
formal Active/Inactive designation. This client uses:

  1. NCBI E-utils esearch to find assays with IC50 data and Active outcomes
     against the gene target (via PubChem BioAssay text search).
  2. PubChem PUG REST to get active CIDs for those assays.
  3. PubChem PUG REST to retrieve compound properties.
"""

from __future__ import annotations

import asyncio
import logging
import os
from typing import Optional

import httpx
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception

from genesis_bio_mcp.models import CompoundActivity, Compounds

logger = logging.getLogger(__name__)

_PUG_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
_NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "genesis-bio-mcp@example.com")

_SEMAPHORE = asyncio.Semaphore(3)


def _is_rate_limited(exc: Exception) -> bool:
    if isinstance(exc, httpx.HTTPStatusError):
        return exc.response.status_code in (429, 503)
    return False


class PubChemClient:
    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client

    async def get_compounds(self, gene_symbol: str) -> Optional[Compounds]:
        """Return active small molecules with bioactivity against a gene target."""
        symbol = gene_symbol.strip().upper()

        # Strategy 1: NCBI Entrez search for BioAssay records with IC50/active data
        aids = await self._search_assays_entrez(symbol)

        # Strategy 2: Fall back to UniProt accession-based AID lookup if Entrez finds nothing
        if not aids:
            aids = await self._get_aids_by_gene_symbol(symbol)

        if not aids:
            logger.info("PubChem: no assays found for gene '%s'", symbol)
            return None

        # Collect active CIDs across assays
        all_compounds: list[CompoundActivity] = []
        for aid in aids[:8]:
            compounds = await self._get_active_cids(aid, symbol)
            all_compounds.extend(compounds)
            if len(all_compounds) >= 100:
                break

        if not all_compounds:
            # Strategy 3: Try fetching compounds that mention the gene via PubChem compound search
            all_compounds = await self._search_compounds_by_name(symbol)

        if not all_compounds:
            return None

        # Deduplicate by CID, keep most potent (lowest IC50)
        by_cid: dict[int, CompoundActivity] = {}
        for c in all_compounds:
            if c.cid not in by_cid:
                by_cid[c.cid] = c
            else:
                existing = by_cid[c.cid]
                if c.activity_value is not None and (
                    existing.activity_value is None or c.activity_value < existing.activity_value
                ):
                    by_cid[c.cid] = c

        active = [c for c in by_cid.values() if c.activity_outcome == "Active"]
        active.sort(key=lambda c: c.activity_value if c.activity_value is not None else float("inf"))

        return Compounds(
            gene_symbol=symbol,
            total_active_compounds=len(active),
            compounds=active[:20],
        )

    async def _search_assays_entrez(self, symbol: str) -> list[int]:
        """Use NCBI Entrez to find BioAssay IDs with active compounds against this gene."""
        term = f'("{symbol}"[Gene Symbol]) AND "active"[Activity Outcome] AND "IC50"[Measurement Type]'
        params = {
            "db": "pcassay",
            "term": term,
            "retmode": "json",
            "retmax": "20",
            "email": _NCBI_EMAIL,
        }
        try:
            # Route through _get() to get retry-on-503 behaviour
            data = await self._get(_ESEARCH_URL, params=params)
            if data is None:
                return []
            ids = data.get("esearchresult", {}).get("idlist", [])
            return [int(i) for i in ids]
        except Exception as exc:
            logger.debug("Entrez BioAssay search failed for '%s': %s", symbol, exc)
            return []

    async def _get_aids_by_gene_symbol(self, symbol: str) -> list[int]:
        """Get assay IDs for a gene symbol from PubChem PUG."""
        url = f"{_PUG_BASE}/assay/target/genesymbol/{symbol}/aids/JSON"
        try:
            data = await self._get(url)
            if data is None:
                return []
            return data.get("IdentifierList", {}).get("AID", [])[:20]
        except Exception as exc:
            logger.debug("PubChem AID lookup failed for '%s': %s", symbol, exc)
            return []

    async def _get_active_cids(self, aid: int, gene_symbol: str) -> list[CompoundActivity]:
        """Get active CIDs for an assay and fetch their properties."""
        url = f"{_PUG_BASE}/assay/aid/{aid}/cids/JSON?cids_type=active"
        try:
            data = await self._get(url)
            if data is None:
                return []
            cids = data.get("IdentifierList", {}).get("CID", [])
            if not cids:
                return []
            return await self._fetch_compound_properties(cids[:20], aid, gene_symbol)
        except Exception as exc:
            logger.debug("PubChem active CIDs failed for AID %d: %s", aid, exc)
            return []

    async def _fetch_compound_properties(
        self, cids: list[int], aid: int, gene_symbol: str
    ) -> list[CompoundActivity]:
        if not cids:
            return []
        cid_str = ",".join(str(c) for c in cids)
        url = f"{_PUG_BASE}/compound/cid/{cid_str}/property/MolecularFormula,MolecularWeight,IUPACName/JSON"
        try:
            data = await self._get(url)
            if data is None:
                return []
            props = data.get("PropertyTable", {}).get("Properties", [])
            return [
                CompoundActivity(
                    cid=p["CID"],
                    name=p.get("IUPACName", f"CID {p['CID']}")[:80],
                    molecular_formula=p.get("MolecularFormula"),
                    molecular_weight=p.get("MolecularWeight"),
                    activity_outcome="Active",
                    assay_id=aid,
                )
                for p in props
                if "CID" in p
            ]
        except Exception as exc:
            logger.debug("PubChem property fetch failed: %s", exc)
            return []

    async def _search_compounds_by_name(self, symbol: str) -> list[CompoundActivity]:
        """Fallback: search for compounds related to the gene by name in PubChem."""
        # Look up well-known drug names via PubChem synonym search
        search_terms = [f"{symbol} inhibitor", symbol]
        for term in search_terms:
            url = f"{_PUG_BASE}/compound/name/{term}/cids/JSON"
            try:
                data = await self._get(url)
                if data is None:
                    continue
                cids = data.get("IdentifierList", {}).get("CID", [])
                if cids:
                    return await self._fetch_compound_properties(cids[:10], aid=0, gene_symbol=symbol)
            except Exception:
                continue
        return []

    @retry(
        retry=retry_if_exception(_is_rate_limited),
        wait=wait_exponential(min=1, max=8),
        stop=stop_after_attempt(3),
    )
    async def _get(self, url: str, params: Optional[dict] = None) -> Optional[dict]:
        async with _SEMAPHORE:
            resp = await self._client.get(url, params=params, timeout=20.0)
            if resp.status_code == 404:
                return None
            resp.raise_for_status()
            return resp.json()
