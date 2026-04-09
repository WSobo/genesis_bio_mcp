"""GWAS Catalog REST API client."""

from __future__ import annotations

import logging
import unicodedata
from typing import Optional

import httpx

from genesis_bio_mcp.models import GwasEvidence, GwasHit

logger = logging.getLogger(__name__)

_BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api"


def _normalize(text: str) -> str:
    """Normalize Unicode to ASCII-comparable form (handles curly apostrophes etc.)."""
    return unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii").lower()


class GwasClient:
    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client

    async def get_evidence(
        self, gene_symbol: str, trait: str, ncbi_gene_id: Optional[str] = None
    ) -> Optional[GwasEvidence]:
        """Return GWAS Catalog hits for a gene–trait pair."""
        symbol = gene_symbol.strip().upper()
        associations = await self._fetch_by_gene_symbol(symbol)

        if not associations and ncbi_gene_id:
            associations = await self._fetch_by_gene_id(ncbi_gene_id)

        if not associations:
            return None

        filtered = _filter_by_trait(associations, trait)
        if not filtered:
            # Return all hits if trait filter yields nothing (gene may have different trait labels)
            filtered = associations

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

    async def _fetch_by_gene_symbol(self, symbol: str) -> list[GwasHit]:
        url = f"{_BASE_URL}/singleNucleotidePolymorphisms/search/findByGene"
        params = {"geneName": symbol, "size": 200}
        return await self._fetch_associations_from_snps(url, params)

    async def _fetch_by_gene_id(self, ncbi_id: str) -> list[GwasHit]:
        url = f"{_BASE_URL}/associations/search/findByEntrezMappedGeneId"
        params = {"entrezMappedGeneId": ncbi_id, "size": 200}
        return await self._fetch_associations(url, params)

    async def _fetch_associations_from_snps(self, url: str, params: dict) -> list[GwasHit]:
        """Fetch SNPs and then resolve to associations."""
        try:
            resp = await self._client.get(url, params=params, timeout=30.0)
            if resp.status_code in (404, 400):
                return []
            resp.raise_for_status()
            body = resp.json()
            snps = body.get("_embedded", {}).get("singleNucleotidePolymorphisms", [])

            hits: list[GwasHit] = []
            for snp in snps[:50]:  # Limit to avoid too many follow-up requests
                assoc_href = snp.get("_links", {}).get("associations", {}).get("href")
                if not assoc_href:
                    continue
                try:
                    ar = await self._client.get(assoc_href, timeout=15.0)
                    ar.raise_for_status()
                    ab = ar.json()
                    for assoc in ab.get("_embedded", {}).get("associations", [])[:3]:
                        hit = _parse_association(assoc)
                        if hit:
                            hits.append(hit)
                except Exception:
                    continue
            return hits
        except Exception as exc:
            logger.warning("GWAS Catalog SNP fetch failed: %s", exc)
            return []

    async def _fetch_associations(self, url: str, params: dict) -> list[GwasHit]:
        try:
            resp = await self._client.get(url, params=params, timeout=30.0)
            if resp.status_code in (404, 400):
                return []
            resp.raise_for_status()
            body = resp.json()
            associations = body.get("_embedded", {}).get("associations", [])
            hits = [_parse_association(a) for a in associations]
            return [h for h in hits if h is not None]
        except Exception as exc:
            logger.warning("GWAS Catalog association fetch failed: %s", exc)
            return []


def _parse_association(assoc: dict) -> Optional[GwasHit]:
    try:
        p_value = assoc.get("pvalue") or assoc.get("pValue")
        if p_value is None:
            return None
        try:
            p_value = float(p_value)
        except (TypeError, ValueError):
            return None

        # Risk allele
        loci = assoc.get("loci", [])
        risk_allele = ""
        mapped_gene = ""
        if loci:
            alleles = loci[0].get("strongestRiskAlleles", [])
            if alleles:
                risk_allele = alleles[0].get("riskAlleleName", "")
            # Mapped genes
            genes = loci[0].get("authorReportedGenes", []) or loci[0].get("entrezMappedGenes", [])
            if genes:
                mapped_gene = genes[0].get("geneName", "")

        # EFO traits
        efo_traits = assoc.get("efoTraits", [])
        trait = efo_traits[0].get("trait", "") if efo_traits else assoc.get("traitName", "")

        # Study info
        study_link = assoc.get("_links", {}).get("study", {}).get("href", "")
        study_accession = ""
        pubmed_id = None
        sample_size = None
        population = None

        # Study may be embedded
        study = assoc.get("study", {})
        if study:
            study_accession = study.get("studyAccession", "")
            pubmed_id = study.get("pubmedId")
            initial_sample = study.get("initialSampleSize", "")
            population = initial_sample[:100] if initial_sample else None

        beta = assoc.get("betaNum") or assoc.get("orPerCopyNum")

        return GwasHit(
            study_accession=study_accession or study_link.split("/")[-1],
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


def _filter_by_trait(hits: list[GwasHit], trait: str) -> list[GwasHit]:
    trait_norm = _normalize(trait)
    return [h for h in hits if trait_norm in _normalize(h.trait)]
