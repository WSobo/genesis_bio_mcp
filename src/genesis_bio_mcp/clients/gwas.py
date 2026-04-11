"""GWAS Catalog REST API client."""

from __future__ import annotations

import asyncio
import logging
import unicodedata

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
        self, gene_symbol: str, trait: str, ncbi_gene_id: str | None = None
    ) -> GwasEvidence | None:
        """Return GWAS Catalog hits for a gene–trait pair."""
        symbol = gene_symbol.strip().upper()

        # Primary: gene-ID path returns fully-embedded study + EFO traits
        associations: list[GwasHit] = []
        if ncbi_gene_id:
            associations = await self._fetch_by_gene_id(ncbi_gene_id)

        # Fallback: SNP path (study links need a follow-up request)
        if not associations:
            associations = await self._fetch_by_gene_symbol(symbol)

        if not associations:
            return None

        # Deduplicate: same association can appear multiple times (multiple SNPs → same locus)
        seen: set[tuple] = set()
        unique: list[GwasHit] = []
        for h in associations:
            key = (h.risk_allele, h.study_accession, h.p_value)
            if key not in seen:
                seen.add(key)
                unique.append(h)
        associations = unique

        filtered = _filter_by_trait(associations, trait)
        if not filtered:
            logger.info(
                "GWAS: no '%s' trait hits for %s — returning None (data gap)",
                trait,
                symbol,
            )
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

    async def _fetch_by_gene_symbol(self, symbol: str) -> list[GwasHit]:
        url = f"{_BASE_URL}/singleNucleotidePolymorphisms/search/findByGene"
        params = {"geneName": symbol, "size": 50}
        return await self._fetch_associations_from_snps(url, params)

    async def _fetch_by_gene_id(self, ncbi_id: str) -> list[GwasHit]:
        url = f"{_BASE_URL}/associations/search/findByEntrezMappedGeneId"
        params = {"entrezMappedGeneId": ncbi_id, "size": 50}
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

            # Fan out all SNP→association fetches concurrently (was sequential)
            results = await asyncio.gather(*[_fetch_snp_assocs(s) for s in snps[:10]])
            raw_assocs = [assoc for batch in results for assoc in batch]

            # Resolve study data for associations that don't embed it
            await self._resolve_study_data(raw_assocs)

            hits = [_parse_association(a) for a in raw_assocs]
            return [h for h in hits if h is not None]
        except Exception as exc:
            logger.warning("GWAS Catalog SNP fetch failed: %s", exc)
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

        # Fan out study fetches concurrently (was sequential)
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
            resp = await self._client.get(url, params=params, timeout=15.0)
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


_TRAIT_SYNONYMS: dict[str, list[str]] = {
    "hypercholesterolemia": ["cholesterol", "ldl", "low-density lipoprotein", "lipid"],
    "hypercholesterolaemia": ["cholesterol", "ldl", "low-density lipoprotein", "lipid"],
    "obesity": [
        "obesity",
        "body mass index",
        "bmi",
        "adiposity",
        "overweight",
        "waist",
    ],
    "rheumatoid arthritis": ["rheumatoid arthritis", "arthritis"],
    "inflammation": ["inflammation", "inflammatory", "c-reactive protein", "crp"],
    "cardiovascular disease": [
        "cardiovascular",
        "coronary artery",
        "myocardial infarction",
        "heart disease",
    ],
    "type 2 diabetes": [
        "type 2 diabetes",
        "t2d",
        "diabetes mellitus",
        "glycated haemoglobin",
        "hba1c",
    ],
    "alzheimer disease": ["alzheimer", "dementia", "cognitive decline"],
    "non-small cell lung carcinoma": [
        "lung cancer",
        "lung carcinoma",
        "non-small cell lung",
    ],
    "squamous cell carcinoma": ["squamous cell", "carcinoma"],
    "pain": ["pain"],
}


def _filter_by_trait(hits: list[GwasHit], trait: str) -> list[GwasHit]:
    trait_norm = _normalize(trait)
    synonyms = _TRAIT_SYNONYMS.get(trait_norm, [])
    match_terms = [trait_norm] + [_normalize(s) for s in synonyms]

    def _matches(hit_trait: str) -> bool:
        ht = _normalize(hit_trait)
        return any(term in ht for term in match_terms)

    return [h for h in hits if _matches(h.trait)]
