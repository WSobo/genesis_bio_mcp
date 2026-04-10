"""Unit tests for API client modules using respx to mock httpx."""

import pytest
import httpx
import respx

from genesis_bio_mcp.clients.uniprot import UniProtClient
from genesis_bio_mcp.clients.open_targets import OpenTargetsClient
from genesis_bio_mcp.clients.depmap import DepMapClient
from genesis_bio_mcp.clients.gwas import GwasClient
from genesis_bio_mcp.clients.pubchem import PubChemClient

from tests.conftest import (
    MOCK_UNIPROT_BRAF,
    MOCK_OT_GENE_SEARCH,
    MOCK_OT_DISEASE_SEARCH,
    MOCK_OT_ASSOCIATION,
    MOCK_DEPMAP_OT_CANCER,
    MOCK_GWAS_ASSOCIATION_RESPONSE,
    MOCK_PUBCHEM_AIDS,
    MOCK_ENTREZ_AIDS,
    MOCK_PUBCHEM_ACTIVE_CIDS,
    MOCK_PUBCHEM_PROPERTIES,
)


# ---------------------------------------------------------------------------
# UniProt client tests
# ---------------------------------------------------------------------------


@respx.mock
async def test_uniprot_get_protein_braf(http_client):
    respx.get("https://rest.uniprot.org/uniprotkb/search").mock(
        return_value=httpx.Response(200, json={"results": [MOCK_UNIPROT_BRAF]})
    )
    client = UniProtClient(http_client)
    result = await client.get_protein("BRAF")

    assert result is not None
    assert result.gene_symbol == "BRAF"
    assert result.uniprot_accession == "P15056"
    assert result.reviewed is True
    assert "kinase" in result.protein_name.lower()
    assert "Cytoplasm" in result.subcellular_locations
    assert any("RAF" in p for p in result.pathways)
    assert len(result.known_variants) == 1
    assert result.known_variants[0].position == "600"
    assert result.known_variants[0].original == "V"
    assert result.known_variants[0].variant == "E"


@respx.mock
async def test_uniprot_returns_none_for_unknown_gene(http_client):
    respx.get("https://rest.uniprot.org/uniprotkb/search").mock(
        return_value=httpx.Response(200, json={"results": []})
    )
    client = UniProtClient(http_client)
    result = await client.get_protein("NOTAREALGENE999")
    assert result is None


@respx.mock
async def test_uniprot_falls_back_to_unreviewed(http_client):
    # First call (reviewed=True) returns empty; second (no filter) returns data
    call_count = 0

    def side_effect(request):
        nonlocal call_count
        call_count += 1
        if call_count == 1:
            return httpx.Response(200, json={"results": []})
        return httpx.Response(200, json={"results": [MOCK_UNIPROT_BRAF]})

    respx.get("https://rest.uniprot.org/uniprotkb/search").mock(side_effect=side_effect)
    client = UniProtClient(http_client)
    result = await client.get_protein("BRAF")
    assert result is not None
    assert call_count == 2


# ---------------------------------------------------------------------------
# Open Targets client tests
# ---------------------------------------------------------------------------


@respx.mock
async def test_open_targets_association_braf_melanoma(http_client):
    # Three sequential GraphQL POST calls: gene search, disease search, association
    route = respx.post("https://api.platform.opentargets.org/api/v4/graphql")
    route.side_effect = [
        httpx.Response(200, json=MOCK_OT_GENE_SEARCH),
        httpx.Response(200, json=MOCK_OT_DISEASE_SEARCH),
        httpx.Response(200, json=MOCK_OT_ASSOCIATION),
    ]

    client = OpenTargetsClient(http_client)
    result = await client.get_association("BRAF", "melanoma")

    assert result is not None
    assert result.gene_symbol == "BRAF"
    assert result.ensembl_id == "ENSG00000157764"
    assert result.disease_efo_id == "EFO_0000389"
    assert result.overall_score == pytest.approx(0.89)
    assert result.somatic_mutation_score == pytest.approx(0.95)
    assert result.known_drug_score == pytest.approx(0.88)
    assert result.evidence_count == 4  # Derived from number of datatypeScores
    assert len(result.evidence_breakdown) == 4


@respx.mock
async def test_open_targets_returns_none_when_gene_not_found(http_client):
    respx.post("https://api.platform.opentargets.org/api/v4/graphql").mock(
        return_value=httpx.Response(200, json={"data": {"search": {"hits": []}}})
    )
    client = OpenTargetsClient(http_client)
    result = await client.get_association("FAKEGENE", "melanoma")
    assert result is None


@respx.mock
async def test_open_targets_returns_none_when_no_association(http_client):
    # Both Bs-filtered query and fallback top-100 return empty rows
    no_assoc = {"data": {"target": {"associatedDiseases": {"count": 0, "rows": []}}}}
    route = respx.post("https://api.platform.opentargets.org/api/v4/graphql")
    route.side_effect = [
        httpx.Response(200, json=MOCK_OT_GENE_SEARCH),
        httpx.Response(200, json=MOCK_OT_DISEASE_SEARCH),
        httpx.Response(200, json=no_assoc),  # Bs filter
        httpx.Response(200, json=no_assoc),  # fallback top-100
    ]
    client = OpenTargetsClient(http_client)
    result = await client.get_association("BRAF", "Alzheimer disease")
    assert result is None


# ---------------------------------------------------------------------------
# DepMap client tests
# ---------------------------------------------------------------------------


@respx.mock
async def test_depmap_returns_essentiality(http_client):
    # DepMap now queries Open Targets cancer associations as proxy
    ot_route = respx.post("https://api.platform.opentargets.org/api/v4/graphql")
    ot_route.side_effect = [
        httpx.Response(200, json=MOCK_OT_GENE_SEARCH),      # gene resolution
        httpx.Response(200, json=MOCK_DEPMAP_OT_CANCER),    # cancer associations
    ]
    client = DepMapClient(http_client, gene_dep_cache={})
    result = await client.get_essentiality("BRAF")

    assert result is not None
    assert result.gene_symbol == "BRAF"
    assert result.mean_ceres_score < 0          # somatic scores mapped to negative proxy
    assert result.fraction_dependent_lines > 0
    assert result.pan_essential is False
    assert len(result.cell_lines) > 0
    assert result.cell_lines[0].is_dependent is True
    assert "Open Targets" in result.data_source


@respx.mock
async def test_depmap_returns_none_when_gene_not_found(http_client):
    respx.post("https://api.platform.opentargets.org/api/v4/graphql").mock(
        return_value=httpx.Response(200, json={"data": {"search": {"hits": []}}})
    )
    client = DepMapClient(http_client, gene_dep_cache={})
    result = await client.get_essentiality("FAKEGENE")
    assert result is None


# ---------------------------------------------------------------------------
# GWAS Catalog client tests
# ---------------------------------------------------------------------------


@respx.mock
async def test_gwas_returns_evidence(http_client):
    # When ncbi_gene_id is provided, findByEntrezMappedGeneId is tried first (primary path).
    # Its response embeds study + efoTraits, so study_accession and trait are populated.
    respx.get(url__regex=r"findByEntrezMappedGeneId").mock(
        return_value=httpx.Response(200, json=MOCK_GWAS_ASSOCIATION_RESPONSE)
    )
    client = GwasClient(http_client)
    result = await client.get_evidence("BRAF", "melanoma", ncbi_gene_id="673")

    assert result is not None
    assert result.gene_symbol == "BRAF"
    assert result.total_associations >= 1
    assert result.associations[0].p_value < 1e-10
    assert result.associations[0].study_accession == "GCST001234"
    assert result.associations[0].trait == "melanoma"


@respx.mock
async def test_gwas_returns_none_when_no_hits(http_client):
    respx.get(url__regex=r"ebi\.ac\.uk/gwas").mock(
        return_value=httpx.Response(200, json={"_embedded": {"singleNucleotidePolymorphisms": [], "associations": []}})
    )
    client = GwasClient(http_client)
    result = await client.get_evidence("FAKEGENE", "nonexistent disease")
    assert result is None


# ---------------------------------------------------------------------------
# PubChem client tests
# ---------------------------------------------------------------------------


@respx.mock
async def test_pubchem_returns_compounds(http_client):
    # New flow: Entrez esearch → active CIDs → compound properties
    respx.get(url__regex=r"esearch\.fcgi").mock(
        return_value=httpx.Response(200, json=MOCK_ENTREZ_AIDS)
    )
    respx.get(url__regex=r"cids/JSON\?cids_type=active").mock(
        return_value=httpx.Response(200, json=MOCK_PUBCHEM_ACTIVE_CIDS)
    )
    respx.get(url__regex=r"property/MolecularFormula").mock(
        return_value=httpx.Response(200, json=MOCK_PUBCHEM_PROPERTIES)
    )
    client = PubChemClient(http_client)
    result = await client.get_compounds("BRAF")

    assert result is not None
    assert result.gene_symbol == "BRAF"
    assert result.total_active_compounds >= 1
    active = [c for c in result.compounds if c.activity_outcome == "Active"]
    assert len(active) >= 1
    assert any("vemurafenib" in c.name.lower() or "dabrafenib" in c.name.lower() for c in active)


@respx.mock
async def test_pubchem_retries_on_503(http_client):
    """Verify tenacity retries on 503 from Entrez."""
    call_count = 0

    def side_effect(request):
        nonlocal call_count
        call_count += 1
        if call_count == 1:
            return httpx.Response(503)
        return httpx.Response(200, json=MOCK_ENTREZ_AIDS)

    respx.get(url__regex=r"esearch\.fcgi").mock(side_effect=side_effect)
    respx.get(url__regex=r"cids/JSON\?cids_type=active").mock(
        return_value=httpx.Response(200, json=MOCK_PUBCHEM_ACTIVE_CIDS)
    )
    respx.get(url__regex=r"property/MolecularFormula").mock(
        return_value=httpx.Response(200, json=MOCK_PUBCHEM_PROPERTIES)
    )
    client = PubChemClient(http_client)
    result = await client.get_compounds("BRAF")

    assert result is not None
    assert call_count == 2  # Retry happened


@respx.mock
async def test_pubchem_returns_none_when_no_assays(http_client):
    # Entrez returns empty, fallback gene symbol AID lookup returns 404,
    # fallback name search also returns 404
    respx.get(url__regex=r"esearch\.fcgi").mock(
        return_value=httpx.Response(200, json={"esearchresult": {"idlist": []}})
    )
    respx.get(url__regex=r"genesymbol").mock(return_value=httpx.Response(404))
    respx.get(url__regex=r"compound/name").mock(return_value=httpx.Response(404))
    client = PubChemClient(http_client)
    result = await client.get_compounds("FAKEGENE")
    assert result is None
