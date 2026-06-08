"""PubChem PUG REST client for compound bioactivity data.

The naive `assay/aid/{aid}/cids/JSON?cids_type=active` endpoint returns 0
CIDs for most BRAF/EGFR/kinase assays because PubChem only flags rows as
"Active" via the per-row Activity Outcome column inside an assay's
*concise* table — the bulk `cids_type=active` index is sparse.

This client therefore:

  1. Resolves the gene symbol to an NCBI GeneID via PubChem's gene domain
     (used to filter panel-assay rows down to the requested target).
  2. Fetches AIDs targeting the gene via
     `assay/target/genesymbol/{symbol}/aids/JSON`.
  3. For each AID, GETs `assay/aid/{aid}/concise/JSON` and parses rows
     where Target GeneID matches and Activity Outcome == "Active".
  4. Dedups by CID, keeps the most potent value, and enriches the top
     results with formula/MW/name via the compound property endpoint.

Activity values in the concise table are reported in µM; we convert to nM
for the model so units match ChEMBL.
"""

from __future__ import annotations

import asyncio
import logging
import os

import httpx
from tenacity import retry, retry_if_exception, stop_after_attempt, wait_exponential

from genesis_bio_mcp.config.settings import settings
from genesis_bio_mcp.models import CompoundActivity, Compounds, SimilarCompound, SimilarCompounds
from genesis_bio_mcp.tools.cheminformatics import tanimoto_similarity

logger = logging.getLogger(__name__)

_PUG_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
_NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "genesis-bio-mcp@example.com")

_SEMAPHORE = asyncio.Semaphore(settings.pubchem_semaphore_limit)

# Concise tables can be ~1 MB each; cap the per-query AID fan-out so
# total payload stays bounded. 10 AIDs against a heavily-drugged gene
# already yields hundreds of unique CIDs after target-row filtering.
_MAX_AIDS_TO_PROBE = 10

# Number of compounds to enrich with formula/MW/IUPAC name. The full
# active set is reported in `total_active_compounds`; only this many are
# materialised in `compounds`.
_MAX_COMPOUNDS_TO_ENRICH = 20

# Human taxonomy ID — PubChem's gene endpoint can return orthologs across
# species; we filter to human to align with the rest of the project.
_HUMAN_TAXONOMY_ID = 9606


def _is_rate_limited(exc: Exception) -> bool:
    if isinstance(exc, httpx.HTTPStatusError):
        return exc.response.status_code in (429, 503)
    return False


class PubChemClient:
    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client

    async def get_compounds(self, gene_symbol: str) -> Compounds | None:
        """Return active small molecules with bioactivity against a gene target."""
        symbol = gene_symbol.strip().upper()

        # ``None`` from an assay-lookup helper signals an UPSTREAM ERROR (timeout /
        # 5xx / network); ``[]`` means PubChem was reachable and genuinely reported
        # no assays. Keeping them apart stops a transient outage from masquerading
        # as a confirmed absence of compounds.
        gene_id, aids = await asyncio.gather(
            self._resolve_gene_id(symbol),
            self._get_aids_by_gene_symbol(symbol),
        )
        primary_errored = aids is None

        if not aids:
            # Primary returned no AIDs (error or genuine empty); try the Entrez fallback.
            fallback = await self._search_assays_entrez(symbol)
            if fallback is None:
                # Fallback also errored. If the primary errored too, every assay-lookup
                # path failed → PubChem is unavailable, not empty.
                if primary_errored:
                    logger.warning("PubChem unavailable for '%s' — assay lookups failed", symbol)
                    return None
                fallback = []  # primary reachably reported no AIDs; trust that
            aids = fallback

        empty = Compounds(gene_symbol=symbol, total_active_compounds=0, compounds=[])
        if not aids:
            logger.info("PubChem: no assays found for gene '%s'", symbol)
            return empty
        if gene_id is None:
            logger.info(
                "PubChem: could not resolve GeneID for '%s' — skipping target filter", symbol
            )

        # Fetch concise tables for the top AIDs concurrently. _get's
        # semaphore enforces the project-wide PubChem rate limit.
        top_aids = aids[:_MAX_AIDS_TO_PROBE]
        results = await asyncio.gather(
            *[self._extract_active_from_concise(aid, gene_id) for aid in top_aids],
            return_exceptions=True,
        )
        all_compounds: list[CompoundActivity] = []
        n_errored = 0
        for r in results:
            if isinstance(r, list):
                all_compounds.extend(r)
            else:
                n_errored += 1

        if not all_compounds:
            # Got AIDs but every concise table failed to load → unavailable, not empty.
            if results and n_errored == len(results):
                logger.warning("PubChem unavailable for '%s' — all assay tables failed", symbol)
                return None
            return empty

        # Dedup by CID, keep the most potent measurement (lowest IC50 in nM).
        by_cid: dict[int, CompoundActivity] = {}
        for c in all_compounds:
            existing = by_cid.get(c.cid)
            if existing is None:
                by_cid[c.cid] = c
                continue
            if c.activity_value is not None and (
                existing.activity_value is None or c.activity_value < existing.activity_value
            ):
                by_cid[c.cid] = c

        active = sorted(
            by_cid.values(),
            key=lambda c: c.activity_value if c.activity_value is not None else float("inf"),
        )

        # Enrich the top compounds with formula/MW/name in a single batch call.
        top = active[:_MAX_COMPOUNDS_TO_ENRICH]
        props = await self._fetch_compound_properties([c.cid for c in top])
        enriched = [
            c.model_copy(
                update={
                    "name": props.get(c.cid, {}).get("name") or c.name,
                    "smiles": props.get(c.cid, {}).get("smiles"),
                    "molecular_formula": props.get(c.cid, {}).get("formula"),
                    "molecular_weight": props.get(c.cid, {}).get("weight"),
                }
            )
            for c in top
        ]

        return Compounds(
            gene_symbol=symbol,
            total_active_compounds=len(active),
            compounds=enriched,
        )

    async def search_similar(
        self,
        smiles: str,
        *,
        mode: str = "similarity",
        threshold: float = 0.9,
        max_results: int = 20,
    ) -> SimilarCompounds | None:
        """Structure-based compound search via PubChem (2D similarity or substructure).

        ``mode="similarity"`` uses PubChem ``fastsimilarity_2d`` (Tanimoto ≥ ``threshold``);
        ``mode="substructure"`` uses ``fastsubstructure``. Hits are enriched with
        name/formula/MW/SMILES, and (similarity mode) ranked by an RDKit Morgan Tanimoto
        computed locally against the query. Returns ``None`` on error.
        """
        query = (smiles or "").strip()
        if not query:
            return None

        if mode == "substructure":
            endpoint = f"{_PUG_BASE}/compound/fastsubstructure/smiles/cids/JSON"
            form = {"smiles": query, "MaxRecords": str(max_results)}
        else:
            mode = "similarity"
            endpoint = f"{_PUG_BASE}/compound/fastsimilarity_2d/smiles/cids/JSON"
            form = {
                "smiles": query,
                "Threshold": str(int(threshold * 100)),
                "MaxRecords": str(max_results),
            }

        # POST the SMILES (path form breaks on '#', '/', etc. in SMILES).
        try:
            async with _SEMAPHORE:
                resp = await self._client.post(endpoint, data=form, timeout=30.0)
            if resp.status_code == 404:
                cids: list[int] = []
            else:
                resp.raise_for_status()
                cids = resp.json().get("IdentifierList", {}).get("CID", []) or []
        except Exception as exc:
            logger.warning("PubChem %s search failed for %s: %s", mode, query, exc)
            return None

        thr = threshold if mode == "similarity" else None
        if not cids:
            return SimilarCompounds(
                query_smiles=query, mode=mode, threshold=thr, total_found=0, hits=[]
            )

        cids = cids[:max_results]
        props = await self._fetch_compound_properties(cids)
        hits: list[SimilarCompound] = []
        for cid in cids:
            p = props.get(cid, {})
            hit_smiles = p.get("smiles")
            tani = (
                tanimoto_similarity(query, hit_smiles)
                if mode == "similarity" and hit_smiles
                else None
            )
            hits.append(
                SimilarCompound(
                    cid=cid,
                    name=p.get("name"),
                    molecular_formula=p.get("formula"),
                    molecular_weight=p.get("weight"),
                    smiles=hit_smiles,
                    tanimoto=tani,
                )
            )
        if mode == "similarity":
            hits.sort(key=lambda h: h.tanimoto if h.tanimoto is not None else -1.0, reverse=True)

        return SimilarCompounds(
            query_smiles=query, mode=mode, threshold=thr, total_found=len(hits), hits=hits
        )

    async def _resolve_gene_id(self, symbol: str) -> int | None:
        """Resolve an HGNC gene symbol to a human NCBI GeneID via PubChem."""
        url = f"{_PUG_BASE}/gene/genesymbol/{symbol}/summary/JSON"
        try:
            data = await self._get(url)
        except Exception as exc:
            logger.debug("PubChem gene resolve failed for '%s': %s", symbol, exc)
            return None
        if not data:
            return None
        summaries = data.get("GeneSummaries", {}).get("GeneSummary", [])
        for gs in summaries:
            if gs.get("TaxonomyID") == _HUMAN_TAXONOMY_ID:
                try:
                    return int(gs["GeneID"])
                except (KeyError, TypeError, ValueError):
                    continue
        # Fall back to the first entry if no explicit human match
        if summaries:
            try:
                return int(summaries[0]["GeneID"])
            except (KeyError, TypeError, ValueError):
                return None
        return None

    async def _get_aids_by_gene_symbol(self, symbol: str) -> list[int] | None:
        """Get all PubChem assay IDs targeting a gene symbol.

        Returns ``[]`` when PubChem is reachable but has no assays for the gene
        (404 / empty list), and ``None`` on an upstream error so the caller can
        tell "no data" apart from "couldn't reach PubChem".
        """
        url = f"{_PUG_BASE}/assay/target/genesymbol/{symbol}/aids/JSON"
        try:
            data = await self._get(url)
            if data is None:
                return []  # 404 → reachable, no assays for this gene
            return data.get("IdentifierList", {}).get("AID", [])
        except Exception as exc:
            logger.debug("PubChem AID lookup failed for '%s': %s", symbol, exc)
            return None  # upstream error — distinct from a genuine empty list

    async def _search_assays_entrez(self, symbol: str) -> list[int] | None:
        """Fallback: NCBI Entrez assay search when PUG REST returns nothing.

        ``None`` on an upstream error, ``[]`` when reachable with no hits.
        """
        term = f'("{symbol}"[Gene Symbol]) AND "active"[Activity Outcome]'
        params = {
            "db": "pcassay",
            "term": term,
            "retmode": "json",
            "retmax": "20",
            "email": _NCBI_EMAIL,
        }
        try:
            data = await self._get(_ESEARCH_URL, params=params)
            if data is None:
                return []
            ids = data.get("esearchresult", {}).get("idlist", [])
            return [int(i) for i in ids]
        except Exception as exc:
            logger.debug("Entrez BioAssay search failed for '%s': %s", symbol, exc)
            return None  # upstream error — distinct from a genuine empty list

    async def _extract_active_from_concise(
        self, aid: int, gene_id: int | None
    ) -> list[CompoundActivity]:
        """Parse an AID's concise bioactivity table; return Active rows for the gene.

        When *gene_id* is None (gene resolution failed), all Active rows are
        returned regardless of target — the caller already restricted the AID
        list by gene symbol so the cross-target leakage is bounded.
        """
        url = f"{_PUG_BASE}/assay/aid/{aid}/concise/JSON"
        # Let _get errors (timeout / 5xx) propagate: the caller gathers with
        # return_exceptions=True and treats "every table errored" as unavailable,
        # not as a genuine zero-compound result. A 404 returns None (no data) below.
        data = await self._get(url)
        if not data:
            return []

        table = data.get("Table", {})
        cols: list[str] = table.get("Columns", {}).get("Column", [])
        rows: list[dict] = table.get("Row", [])
        if not cols or not rows:
            return []

        try:
            cid_idx = cols.index("CID")
            outcome_idx = cols.index("Activity Outcome")
            value_idx = cols.index("Activity Value [uM]")
            gene_idx = cols.index("Target GeneID") if "Target GeneID" in cols else None
            name_idx = cols.index("Activity Name") if "Activity Name" in cols else None
        except ValueError:
            return []

        gene_id_str = str(gene_id) if gene_id is not None else None
        out: list[CompoundActivity] = []
        for row in rows:
            cells = row.get("Cell", [])
            if len(cells) <= max(cid_idx, outcome_idx, value_idx):
                continue
            if cells[outcome_idx] != "Active":
                continue
            if gene_idx is not None and gene_id_str is not None:
                if cells[gene_idx] != gene_id_str:
                    continue
            cid_raw = cells[cid_idx]
            if not cid_raw:
                continue
            try:
                cid = int(cid_raw)
            except (TypeError, ValueError):
                continue
            value_um = (cells[value_idx] or "").strip()
            try:
                value_nm = float(value_um) * 1000.0 if value_um else None
            except ValueError:
                value_nm = None
            activity_type: str | None = None
            if name_idx is not None and name_idx < len(cells):
                activity_type = (cells[name_idx] or "").strip() or None
            out.append(
                CompoundActivity(
                    cid=cid,
                    name=f"CID {cid}",
                    activity_outcome="Active",
                    activity_value=value_nm,
                    activity_type=activity_type,
                    assay_id=aid,
                )
            )
        return out

    async def _fetch_compound_properties(self, cids: list[int]) -> dict[int, dict]:
        """Batch-fetch compound properties; return ``{cid: {formula, weight, name}}``."""
        if not cids:
            return {}
        cid_str = ",".join(str(c) for c in cids)
        # NOTE: PubChem deprecated `CanonicalSMILES`; the current property is
        # `SMILES` (it also still returns `ConnectivitySMILES`). Confirmed
        # against the live PUG REST schema.
        url = (
            f"{_PUG_BASE}/compound/cid/{cid_str}"
            "/property/MolecularFormula,MolecularWeight,IUPACName,SMILES/JSON"
        )
        try:
            data = await self._get(url)
        except Exception as exc:
            logger.debug("PubChem property batch fetch failed: %s", exc)
            return {}
        if not data:
            return {}
        props = data.get("PropertyTable", {}).get("Properties", [])
        out: dict[int, dict] = {}
        for p in props:
            if "CID" not in p:
                continue
            try:
                cid = int(p["CID"])
            except (TypeError, ValueError):
                continue
            mw_raw = p.get("MolecularWeight")
            try:
                mw = float(mw_raw) if mw_raw is not None else None
            except (TypeError, ValueError):
                mw = None
            name = p.get("IUPACName")
            # PubChem returns the structure under `SMILES`; older payloads used
            # `CanonicalSMILES`. Accept either so the parser is release-robust.
            smiles = p.get("SMILES") or p.get("CanonicalSMILES")
            out[cid] = {
                "formula": p.get("MolecularFormula"),
                "weight": mw,
                "name": (name[:80] if isinstance(name, str) else None),
                "smiles": smiles if isinstance(smiles, str) else None,
            }
        return out

    @retry(
        retry=retry_if_exception(_is_rate_limited),
        wait=wait_exponential(min=1, max=8),
        stop=stop_after_attempt(3),
    )
    async def _get(self, url: str, params: dict | None = None) -> dict | None:
        async with _SEMAPHORE:
            resp = await self._client.get(url, params=params, timeout=20.0)
            if resp.status_code == 404:
                return None
            resp.raise_for_status()
            return resp.json()
