"""Regression tests: GWAS uses the cheap gene-ID path first and only falls back
to the expensive SNP fan-out when the gene-ID path is empty.

Background: in a multi-gene compare_targets burst, every gene gapped out of GWAS
even though each returned data when queried alone. Cause: the SNP path
(findByGene -> per-SNP associations -> study details) fires up to ~36
semaphore-gated calls per gene; running it for all genes saturated the shared
GWAS semaphore and starved each gene's single cheap gene-ID call. Fix: make the
SNP fan-out a fallback so the common case is one call per gene.
"""

from __future__ import annotations

import httpx
import respx

from genesis_bio_mcp.clients.gwas import GwasClient

_GENE_ID_HITS = {
    "_embedded": {
        "associations": [
            {
                "pvalue": 4e-20,
                "loci": [
                    {
                        "strongestRiskAlleles": [{"riskAlleleName": "rs11591147-T"}],
                        "authorReportedGenes": [{"geneName": "PCSK9"}],
                    }
                ],
                "efoTraits": [{"trait": "LDL cholesterol levels", "uri": None}],
            }
        ]
    }
}
_EMPTY_GENE_ID = {"_embedded": {"associations": []}}
_EMPTY_SNP = {"_embedded": {"singleNucleotidePolymorphisms": []}}


@respx.mock
async def test_gene_id_hits_skip_the_snp_fanout(monkeypatch, tmp_path) -> None:
    """When the gene-ID path returns hits, the expensive SNP fan-out is not called."""
    monkeypatch.setattr(
        "genesis_bio_mcp.config.settings.settings.gwas_cache_path", tmp_path / "gwas.json"
    )
    gene_id_route = respx.get(url__regex=r"findByEntrezMappedGeneId").mock(
        return_value=httpx.Response(200, json=_GENE_ID_HITS)
    )
    snp_route = respx.get(url__regex=r"findByGene").mock(
        return_value=httpx.Response(200, json=_EMPTY_SNP)
    )

    async with httpx.AsyncClient() as http:
        client = GwasClient(http, efo_resolver=None)
        result = await client.get_evidence("PCSK9", "hypercholesterolemia", ncbi_gene_id="255738")

    assert result is not None  # gene-level hits surfaced
    assert gene_id_route.called
    assert snp_route.call_count == 0  # the ~36-call fan-out was skipped


@respx.mock
async def test_empty_gene_id_falls_back_to_snp(monkeypatch, tmp_path) -> None:
    """When the gene-ID path is empty, the SNP fan-out still runs as a fallback."""
    monkeypatch.setattr(
        "genesis_bio_mcp.config.settings.settings.gwas_cache_path", tmp_path / "gwas.json"
    )
    gene_id_route = respx.get(url__regex=r"findByEntrezMappedGeneId").mock(
        return_value=httpx.Response(200, json=_EMPTY_GENE_ID)
    )
    snp_route = respx.get(url__regex=r"findByGene").mock(
        return_value=httpx.Response(200, json=_EMPTY_SNP)
    )

    async with httpx.AsyncClient() as http:
        client = GwasClient(http, efo_resolver=None)
        result = await client.get_evidence("FOO", "bar", ncbi_gene_id="999")

    assert gene_id_route.called
    assert snp_route.called  # fallback engaged because the gene-ID path was empty
    assert result is None  # both empty -> genuine data gap


@respx.mock
async def test_no_gene_id_uses_snp_path(monkeypatch, tmp_path) -> None:
    """With no NCBI gene ID, the SNP path is the only option."""
    monkeypatch.setattr(
        "genesis_bio_mcp.config.settings.settings.gwas_cache_path", tmp_path / "gwas.json"
    )
    gene_id_route = respx.get(url__regex=r"findByEntrezMappedGeneId").mock(
        return_value=httpx.Response(200, json=_EMPTY_GENE_ID)
    )
    snp_route = respx.get(url__regex=r"findByGene").mock(
        return_value=httpx.Response(200, json=_EMPTY_SNP)
    )

    async with httpx.AsyncClient() as http:
        client = GwasClient(http, efo_resolver=None)
        await client.get_evidence("FOO", "bar")  # no ncbi_gene_id

    assert gene_id_route.call_count == 0  # gene-ID path skipped without an id
    assert snp_route.called
