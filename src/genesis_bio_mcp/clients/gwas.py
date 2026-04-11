"""GWAS Catalog REST API client."""

from __future__ import annotations

import asyncio
import json
import logging
import time
from pathlib import Path

import httpx

from genesis_bio_mcp.config.trait_synonyms import filter_by_trait
from genesis_bio_mcp.models import GwasEvidence, GwasHit

logger = logging.getLogger(__name__)

_BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api"
_CACHE_PATH = Path("data/gwas_cache.json")
_CACHE_TTL_SECS = 86400  # 24 hours


def _load_gwas_cache() -> dict[str, dict]:
    try:
        if _CACHE_PATH.exists():
            return json.loads(_CACHE_PATH.read_text())
    except Exception as exc:
        logger.warning("Failed to load GWAS cache: %s", repr(exc))
    return {}


def _get_cached(cache: dict[str, dict], symbol: str, trait: str) -> GwasEvidence | None:
    entry = cache.get(f"{symbol}:{trait}")
    if not entry:
        return None
    if time.time() - entry.get("fetched_at", 0) > _CACHE_TTL_SECS:
        return None
    try:
        return GwasEvidence.model_validate(entry["result"])
    except Exception:
        return None


def _set_cached(cache: dict[str, dict], symbol: str, trait: str, result: GwasEvidence) -> None:
    cache[f"{symbol}:{trait}"] = {"result": result.model_dump(), "fetched_at": time.time()}
    try:
        _CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
        _CACHE_PATH.write_text(json.dumps(cache, indent=2, default=str))
    except Exception as exc:
        logger.warning("Failed to write GWAS cache: %s", repr(exc))


def _process_for_trait(raw: list[GwasHit], symbol: str, trait: str) -> GwasEvidence | None:
    """Dedup, trait-filter, and package raw hits into a GwasEvidence result."""
    # Deduplicate: same association can appear via both fetch paths
    seen: set[tuple] = set()
    unique: list[GwasHit] = []
    for h in raw:
        key = (h.risk_allele, h.study_accession, h.p_value)
        if key not in seen:
            seen.add(key)
            unique.append(h)

    filtered = filter_by_trait(unique, trait)
    if not filtered:
        return None

    # Remove zero p-value artifacts and sort
    filtered = [a for a in filtered if a.p_value and a.p_value > 0]
    filtered.sort(key=lambda h: h.p_value)

    if not filtered:
        return None

    return GwasEvidence(
        gene_symbol=symbol,
        trait_query=trait,
        total_associations=len(filtered),
        associations=filtered[:50],
        strongest_p_value=filtered[0].p_value if filtered else None,
    )


class GwasClient:
    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client
        self._disk_cache: dict[str, dict] = _load_gwas_cache()
        # Session-level gene cache: raw associations per gene symbol.
        # Avoids re-fetching the same gene for different trait queries
        # (e.g. PTGS2 for "inflammation" then COX2→PTGS2 for "pain").
        self._gene_cache: dict[str, list[GwasHit]] = {}

    async def get_evidence(
        self, gene_symbol: str, trait: str, ncbi_gene_id: str | None = None
    ) -> GwasEvidence | None:
        """Return GWAS Catalog hits for a gene–trait pair."""
        symbol = gene_symbol.strip().upper()

        # Session cache: fetch raw associations once per gene, filter per trait
        if symbol not in self._gene_cache:
            self._gene_cache[symbol] = await self._fetch_all(symbol, ncbi_gene_id)

        raw = self._gene_cache[symbol]

        if not raw:
            # Disk cache fallback — timeouts are infrastructure failures,
            # not evidence of absent GWAS signal.
            cached = _get_cached(self._disk_cache, symbol, trait)
            if cached is not None:
                logger.info("GWAS live fetch empty for %s/%s, serving cached result", symbol, trait)
                return cached
            return None

        result = _process_for_trait(raw, symbol, trait)
        if result:
            _set_cached(self._disk_cache, symbol, trait, result)
        else:
            logger.info(
                "GWAS: no '%s' trait hits for %s — returning None (data gap)",
                trait,
                symbol,
            )
        return result

    async def _fetch_all(self, symbol: str, ncbi_gene_id: str | None) -> list[GwasHit]:
        """Run primary (gene-ID) and SNP paths concurrently, return merged results.

        Uses asyncio.wait with a 15s global bound so worst case = max(paths),
        not sum(paths). Whatever completes within the window is kept; the rest
        is cancelled.
        """
        tasks: list[asyncio.Task] = []
        if ncbi_gene_id:
            tasks.append(asyncio.create_task(self._fetch_by_gene_id(ncbi_gene_id)))
        tasks.append(asyncio.create_task(self._fetch_by_gene_symbol(symbol)))

        if not tasks:
            return []

        # 15s global bound: if primary times out at 25s but SNP returns in 12s,
        # we get SNP results at the 15s mark instead of waiting 45s total.
        done, pending = await asyncio.wait(tasks, timeout=15.0)
        for t in pending:
            t.cancel()

        merged: list[GwasHit] = []
        for t in done:
            try:
                merged.extend(t.result())
            except Exception as exc:
                logger.warning("GWAS fetch path failed: %s", repr(exc))

        if not merged and pending:
            logger.warning("GWAS: all fetch paths timed out for %s within 15s window", symbol)

        return merged

    async def _fetch_by_gene_symbol(self, symbol: str) -> list[GwasHit]:
        url = f"{_BASE_URL}/singleNucleotidePolymorphisms/search/findByGene"
        params = {"geneName": symbol, "size": 50}
        return await self._fetch_associations_from_snps(url, params)

    async def _fetch_by_gene_id(self, ncbi_id: str) -> list[GwasHit]:
        url = f"{_BASE_URL}/associations/search/findByEntrezMappedGeneId"
        # size=200: primary path is a single HTTP call regardless of result count;
        # cutting this caused PCSK9 to lose 12 real lipid-trait hits.
        params = {"entrezMappedGeneId": ncbi_id, "size": 200}
        return await self._fetch_associations(url, params)

    async def _fetch_associations_from_snps(self, url: str, params: dict) -> list[GwasHit]:
        """Fetch SNPs, follow their association links concurrently, resolve study data."""
        try:
            resp = await self._client.get(url, params=params, timeout=15.0)
            if resp.status_code in (404, 400):
                return []
            resp.raise_for_status()
            body = resp.json()
            snps = body.get("_embedded", {}).get("singleNucleotidePolymorphisms", [])

            async def _fetch_snp_assocs(snp: dict) -> list[dict]:
                assoc_href = snp.get("_links", {}).get("associations", {}).get("href")
                if not assoc_href:
                    return []
                try:
                    ar = await self._client.get(assoc_href, timeout=8.0)
                    ar.raise_for_status()
                    ab = ar.json()
                    return ab.get("_embedded", {}).get("associations", [])[:3]
                except Exception:
                    return []

            # Fan out all SNP→association fetches concurrently
            results = await asyncio.gather(*[_fetch_snp_assocs(s) for s in snps[:10]])
            raw_assocs = [assoc for batch in results for assoc in batch]

            # Resolve study data for associations that don't embed it
            await self._resolve_study_data(raw_assocs)

            hits = [_parse_association(a) for a in raw_assocs]
            return [h for h in hits if h is not None]
        except Exception as exc:
            logger.warning("GWAS Catalog SNP fetch failed: %s", repr(exc))
            return []

    async def _resolve_study_data(self, assocs: list[dict]) -> None:
        """Follow study sub-resource links for associations missing embedded study data.

        GWAS HAL associations from the SNP path don't embed the study object; the link
        is ``associations/{id}/study`` (not ``studies/GCST...``). We batch-fetch unique
        study URLs and inject the result back into the raw assoc dicts so _parse_association
        can extract studyAccession and diseaseTrait.
        """
        # Collect unique study URLs that need resolution
        pending: dict[str, None] = {}  # ordered set
        for assoc in assocs:
            if assoc.get("study"):
                continue
            study_link = assoc.get("_links", {}).get("study", {}).get("href", "")
            if study_link:
                pending[study_link] = None

        async def _fetch_study(link: str) -> tuple[str, dict]:
            try:
                sr = await self._client.get(link, timeout=5.0)
                sr.raise_for_status()
                return link, sr.json()
            except Exception:
                return link, {}

        # Fan out study fetches concurrently
        fetched = await asyncio.gather(*[_fetch_study(link) for link in pending])
        study_cache = dict(fetched)

        for assoc in assocs:
            if assoc.get("study"):
                continue
            study_link = assoc.get("_links", {}).get("study", {}).get("href", "")
            if study_link and study_link in study_cache:
                assoc["study"] = study_cache[study_link]

    async def _fetch_associations(self, url: str, params: dict) -> list[GwasHit]:
        try:
            # 25s: primary path is one HTTP call; GWAS Catalog can be slow for
            # gene IDs with many associations (TNF, FTO were timing out at 15s).
            resp = await self._client.get(url, params=params, timeout=25.0)
            if resp.status_code in (404, 400):
                return []
            resp.raise_for_status()
            body = resp.json()
            associations = body.get("_embedded", {}).get("associations", [])
            hits = [_parse_association(a) for a in associations]
            return [h for h in hits if h is not None]
        except Exception as exc:
            logger.warning("GWAS Catalog association fetch failed: %s", repr(exc))
            return []


def _parse_association(assoc: dict) -> GwasHit | None:
    try:
        p_value = assoc.get("pvalue") or assoc.get("pValue")
        if p_value is None:
            return None
        try:
            p_value = float(p_value)
        except (TypeError, ValueError):
            return None

        # Risk allele and mapped gene
        loci = assoc.get("loci", [])
        risk_allele = ""
        mapped_gene = ""
        if loci:
            alleles = loci[0].get("strongestRiskAlleles", [])
            if alleles:
                risk_allele = alleles[0].get("riskAlleleName", "")
            genes = loci[0].get("authorReportedGenes", []) or loci[0].get("entrezMappedGenes", [])
            if genes:
                mapped_gene = genes[0].get("geneName", "")

        # EFO traits — present in gene-ID-path responses; absent in SNP-path responses
        # (SNP-path study objects embed diseaseTrait instead)
        efo_traits = assoc.get("efoTraits", [])
        if efo_traits:
            trait = efo_traits[0].get("trait", "")
        else:
            # Fall back to study.diseaseTrait.trait (populated after _resolve_study_data)
            trait = assoc.get("study", {}).get("diseaseTrait", {}).get("trait", "") or assoc.get(
                "traitName", ""
            )

        # Study info — may be embedded directly or injected by _resolve_study_data
        study = assoc.get("study", {})
        study_accession = study.get("studyAccession", "")
        pubmed_id = study.get("pubmedId")
        initial_sample = study.get("initialSampleSize", "")
        population = initial_sample[:100] if initial_sample else None
        sample_size = None

        # Final fallback: extract GCST from study link if link points directly to studies/GCST*
        if not study_accession:
            study_link = assoc.get("_links", {}).get("study", {}).get("href", "")
            parts = study_link.rstrip("/").split("/")
            candidate = parts[-1] if parts else ""
            if candidate.upper().startswith("GCST"):
                study_accession = candidate

        beta = assoc.get("betaNum") or assoc.get("orPerCopyNum")

        return GwasHit(
            study_accession=study_accession,
            trait=trait,
            mapped_gene=mapped_gene,
            risk_allele=risk_allele,
            p_value=p_value,
            beta_or_or=float(beta) if beta is not None else None,
            sample_size=sample_size,
            population=population,
            pubmed_id=str(pubmed_id) if pubmed_id else None,
        )
    except Exception as exc:
        logger.debug("Failed to parse GWAS association: %s", exc)
        return None
