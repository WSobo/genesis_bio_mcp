"""Unit tests for API client modules using respx to mock httpx."""

import httpx
import pytest
import respx

from genesis_bio_mcp.clients.alphafold import AlphaFoldClient
from genesis_bio_mcp.clients.biogrid import BioGRIDClient
from genesis_bio_mcp.clients.clinical_trials import ClinicalTrialsClient
from genesis_bio_mcp.clients.depmap import DepMapClient
from genesis_bio_mcp.clients.dgidb import DGIdbClient
from genesis_bio_mcp.clients.gnomad import GnomADClient
from genesis_bio_mcp.clients.gwas import GwasClient
from genesis_bio_mcp.clients.iedb import IEDBClient
from genesis_bio_mcp.clients.interpro import InterProClient
from genesis_bio_mcp.clients.mavedb import MaveDBClient
from genesis_bio_mcp.clients.myvariant import MyVariantClient
from genesis_bio_mcp.clients.open_targets import OpenTargetsClient
from genesis_bio_mcp.clients.pubchem import PubChemClient
from genesis_bio_mcp.clients.reactome import ReactomeClient
from genesis_bio_mcp.clients.sabdab import SAbDabClient
from genesis_bio_mcp.clients.string_db import StringDbClient
from genesis_bio_mcp.clients.uniprot import UniProtClient
from genesis_bio_mcp.config.efo_resolver import EFOResolver, EFOTerm
from genesis_bio_mcp.config.trait_synonyms import filter_by_trait
from genesis_bio_mcp.models import GwasHit
from tests.conftest import (
    MOCK_DEPMAP_OT_CANCER,
    MOCK_ENTREZ_AIDS,
    MOCK_GWAS_ASSOCIATION_RESPONSE,
    MOCK_OT_ASSOCIATION,
    MOCK_OT_DISEASE_SEARCH,
    MOCK_OT_GENE_SEARCH,
    MOCK_PUBCHEM_ACTIVE_CIDS,
    MOCK_PUBCHEM_PROPERTIES,
    MOCK_UNIPROT_BRAF,
    MOCK_UNIPROT_FASTA_BRAF,
    MOCK_UNIPROT_FASTA_BRAF_SEQUENCE,
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
async def test_uniprot_parses_disulfide_bonds(http_client):
    respx.get("https://rest.uniprot.org/uniprotkb/search").mock(
        return_value=httpx.Response(200, json={"results": [MOCK_UNIPROT_BRAF]})
    )
    client = UniProtClient(http_client)
    result = await client.get_protein("BRAF")
    assert result is not None
    # BRAF mock has one DISULFID feature at positions 157 and 162
    assert result.disulfide_bond_positions == [157, 162]


@respx.mock
async def test_uniprot_get_sequence_happy_path(http_client):
    respx.get("https://rest.uniprot.org/uniprotkb/P15056.fasta").mock(
        return_value=httpx.Response(200, text=MOCK_UNIPROT_FASTA_BRAF)
    )
    client = UniProtClient(http_client)
    result = await client.get_sequence("P15056")
    assert result is not None
    seq, organism, description = result
    assert seq == MOCK_UNIPROT_FASTA_BRAF_SEQUENCE
    assert "Homo sapiens" in organism
    assert "B-raf" in description


@respx.mock
async def test_uniprot_get_sequence_slice(http_client):
    respx.get("https://rest.uniprot.org/uniprotkb/P15056.fasta").mock(
        return_value=httpx.Response(200, text=MOCK_UNIPROT_FASTA_BRAF)
    )
    client = UniProtClient(http_client)
    result = await client.get_sequence("P15056", start=1, end=10)
    assert result is not None
    seq, _, _ = result
    assert seq == MOCK_UNIPROT_FASTA_BRAF_SEQUENCE[:10]


@respx.mock
async def test_uniprot_get_sequence_404(http_client):
    respx.get("https://rest.uniprot.org/uniprotkb/BADACC.fasta").mock(
        return_value=httpx.Response(404, text="Not found")
    )
    client = UniProtClient(http_client)
    result = await client.get_sequence("BADACC")
    assert result is None


@respx.mock
async def test_uniprot_get_sequence_network_error(http_client):
    respx.get("https://rest.uniprot.org/uniprotkb/P15056.fasta").mock(
        side_effect=httpx.ConnectError("boom")
    )
    client = UniProtClient(http_client)
    result = await client.get_sequence("P15056")
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
        httpx.Response(200, json=MOCK_OT_GENE_SEARCH),  # gene resolution
        httpx.Response(200, json=MOCK_DEPMAP_OT_CANCER),  # cancer associations
    ]
    client = DepMapClient(http_client, gene_dep_cache={})
    result = await client.get_essentiality("BRAF")

    assert result is not None
    assert result.gene_symbol == "BRAF"
    assert result.mean_ceres_score < 0  # somatic scores mapped to negative proxy
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
        return_value=httpx.Response(
            200,
            json={"_embedded": {"singleNucleotidePolymorphisms": [], "associations": []}},
        )
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


# ---------------------------------------------------------------------------
# AlphaFold + RCSB PDB client tests
# ---------------------------------------------------------------------------


_MOCK_ALPHAFOLD_RESPONSE = [
    {
        "uniprotAccession": "P15056",
        "meanPlddt": 91.5,
        "pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P15056-F1-model_v4.pdb",
        "latestVersion": 4,
    }
]

_MOCK_RCSB_SEARCH = {
    "total_count": 2,
    "result_set": [
        {"identifier": "4MNE"},
        {"identifier": "3OG7"},
    ],
}

_MOCK_RCSB_ENTRY = {
    "exptl": [{"method": "X-RAY DIFFRACTION"}],
    "refine": [{"ls_d_res_high": 2.1}],
    "rcsb_entry_info": {"nonpolymer_entity_count": 1},
    "rcsb_accession_info": {"deposit_date": "2012-03-15"},
}


@respx.mock
async def test_alphafold_get_structure(http_client):
    respx.get(url__regex=r"alphafold\.ebi\.ac\.uk/api/prediction").mock(
        return_value=httpx.Response(200, json=_MOCK_ALPHAFOLD_RESPONSE)
    )
    respx.post(url__regex=r"rcsbsearch").mock(
        return_value=httpx.Response(200, json=_MOCK_RCSB_SEARCH)
    )
    respx.get(url__regex=r"data\.rcsb\.org/rest/v1/core/entry").mock(
        return_value=httpx.Response(200, json=_MOCK_RCSB_ENTRY)
    )
    client = AlphaFoldClient(http_client)
    result = await client.get_structure("BRAF", uniprot_accession="P15056")

    assert result is not None
    assert result.gene_symbol == "BRAF"
    assert result.alphafold_plddt == pytest.approx(91.5)
    assert result.alphafold_version == "v4"
    assert result.total_pdb_structures == 2
    assert result.has_ligand_bound is True
    assert result.best_resolution == pytest.approx(2.1)
    assert len(result.experimental_structures) == 2
    assert result.experimental_structures[0].pdb_id in ("4MNE", "3OG7")


@respx.mock
async def test_alphafold_returns_none_without_accession(http_client):
    client = AlphaFoldClient(http_client)
    result = await client.get_structure("BRAF", uniprot_accession=None)
    assert result is None


@respx.mock
async def test_alphafold_handles_404(http_client):
    respx.get(url__regex=r"alphafold\.ebi\.ac\.uk/api/prediction").mock(
        return_value=httpx.Response(404)
    )
    respx.post(url__regex=r"rcsbsearch").mock(return_value=httpx.Response(204))
    client = AlphaFoldClient(http_client)
    result = await client.get_structure("UNKNOWNGENE", uniprot_accession="Q99999")
    assert result is not None
    assert result.alphafold_plddt is None
    assert result.total_pdb_structures == 0


# ---------------------------------------------------------------------------
# STRING client tests
# ---------------------------------------------------------------------------


_MOCK_STRING_RESOLVE = [{"stringId": "9606.ENSP00000288602", "preferredName": "BRAF"}]

_MOCK_STRING_NETWORK = [
    {
        "stringId_A": "9606.ENSP00000288602",
        "preferredName_A": "BRAF",
        "stringId_B": "9606.ENSP00000395687",
        "preferredName_B": "MAP2K1",
        "score": 997,
        "escore": 400,
        "dscore": 900,
        "cscore": 0,
        "tscore": 0,
        "hscore": 0,
        "ascore": 0,
        "fscore": 0,
    },
    {
        "stringId_A": "9606.ENSP00000288602",
        "preferredName_A": "BRAF",
        "stringId_B": "9606.ENSP00000410194",
        "preferredName_B": "RAF1",
        "score": 962,
        "escore": 300,
        "dscore": 800,
        "cscore": 0,
        "tscore": 0,
        "hscore": 0,
        "ascore": 0,
        "fscore": 0,
    },
]


@respx.mock
async def test_string_get_interactome(http_client):
    respx.get(url__regex=r"string-db\.org/api/json/get_string_ids").mock(
        return_value=httpx.Response(200, json=_MOCK_STRING_RESOLVE)
    )
    respx.get(url__regex=r"string-db\.org/api/json/network").mock(
        return_value=httpx.Response(200, json=_MOCK_STRING_NETWORK)
    )
    client = StringDbClient(http_client)
    result = await client.get_interactome("BRAF")

    assert result is not None
    assert result.gene_symbol == "BRAF"
    assert result.total_partners == 2
    symbols = [i.gene_symbol for i in result.top_interactors]
    assert "MAP2K1" in symbols
    assert "RAF1" in symbols
    # Sorted descending by score
    assert result.top_interactors[0].score >= result.top_interactors[1].score


@respx.mock
async def test_string_returns_none_when_unresolvable(http_client):
    respx.get(url__regex=r"string-db\.org/api/json/get_string_ids").mock(
        return_value=httpx.Response(200, json=[])
    )
    client = StringDbClient(http_client)
    result = await client.get_interactome("FAKEGENE")
    assert result is None


@respx.mock
async def test_string_returns_empty_interactors_on_network_failure(http_client):
    """Resolve succeeds but network fetch fails — returns ProteinInteractome with empty list."""
    respx.get(url__regex=r"string-db\.org/api/json/get_string_ids").mock(
        return_value=httpx.Response(200, json=_MOCK_STRING_RESOLVE)
    )
    respx.get(url__regex=r"string-db\.org/api/json/network").mock(return_value=httpx.Response(500))
    client = StringDbClient(http_client)
    result = await client.get_interactome("BRAF")

    assert result is not None
    assert result.gene_symbol == "BRAF"
    assert result.top_interactors == []
    assert result.total_partners == 0


# ---------------------------------------------------------------------------
# BioGRID client tests
# ---------------------------------------------------------------------------

_MOCK_BIOGRID_RESPONSE = {
    "12345": {
        "OFFICIAL_SYMBOL_A": "BRAF",
        "OFFICIAL_SYMBOL_B": "RAF1",
        "EXPERIMENTAL_SYSTEM": "Two-hybrid",
        "EXPERIMENTAL_SYSTEM_TYPE": "physical",
        "PUBMED_ID": "12345678",
        "THROUGHPUT": "Low Throughput",
    },
    "12346": {
        "OFFICIAL_SYMBOL_A": "BRAF",
        "OFFICIAL_SYMBOL_B": "MAP2K1",
        "EXPERIMENTAL_SYSTEM": "Co-immunoprecipitation",
        "EXPERIMENTAL_SYSTEM_TYPE": "physical",
        "PUBMED_ID": "23456789",
        "THROUGHPUT": "Low Throughput",
    },
    "12347": {
        "OFFICIAL_SYMBOL_A": "RAF1",
        "OFFICIAL_SYMBOL_B": "BRAF",
        "EXPERIMENTAL_SYSTEM": "Two-hybrid",
        "EXPERIMENTAL_SYSTEM_TYPE": "physical",
        "PUBMED_ID": "34567890",
        "THROUGHPUT": "Low Throughput",
    },
}


@respx.mock
async def test_biogrid_get_interactions_happy_path(http_client, monkeypatch):
    monkeypatch.setenv("BIOGRID_ACCESS_KEY", "test-key")
    respx.get(url__regex=r"webservice\.thebiogrid\.org/interactions").mock(
        return_value=httpx.Response(200, json=_MOCK_BIOGRID_RESPONSE)
    )
    client = BioGRIDClient(http_client)
    result = await client.get_interactions("BRAF")

    assert result is not None
    assert result.gene_symbol == "BRAF"
    assert result.total_interactions == 3
    assert result.unique_partners >= 2
    partners = {ix.interactor_a for ix in result.interactions} | {
        ix.interactor_b for ix in result.interactions
    }
    assert "RAF1" in partners
    assert "MAP2K1" in partners


@respx.mock
async def test_biogrid_returns_none_without_api_key(http_client, monkeypatch):
    monkeypatch.delenv("BIOGRID_ACCESS_KEY", raising=False)
    client = BioGRIDClient(http_client)
    result = await client.get_interactions("BRAF")
    assert result is None


@respx.mock
async def test_biogrid_returns_empty_on_403(http_client, monkeypatch):
    monkeypatch.setenv("BIOGRID_ACCESS_KEY", "bad-key")
    respx.get(url__regex=r"webservice\.thebiogrid\.org/interactions").mock(
        return_value=httpx.Response(403, json={"error": "Forbidden"})
    )
    client = BioGRIDClient(http_client)
    result = await client.get_interactions("BRAF")
    assert result is None


@respx.mock
async def test_biogrid_empty_results(http_client, monkeypatch):
    monkeypatch.setenv("BIOGRID_ACCESS_KEY", "test-key")
    respx.get(url__regex=r"webservice\.thebiogrid\.org/interactions").mock(
        return_value=httpx.Response(200, json={})
    )
    client = BioGRIDClient(http_client)
    result = await client.get_interactions("UNKNOWNGENE")

    assert result is not None
    assert result.total_interactions == 0
    assert result.unique_partners == 0
    assert result.interactions == []


# ---------------------------------------------------------------------------
# DGIdb client tests
# ---------------------------------------------------------------------------


_MOCK_DGIDB_RESPONSE = {
    "data": {
        "genes": {
            "nodes": [
                {
                    "name": "BRAF",
                    "interactions": [
                        {
                            "drug": {"name": "Vemurafenib", "approved": True},
                            "interactionTypes": [
                                {"type": "inhibitor", "directionality": "inhibitory"}
                            ],
                            "interactionClaims": [{"source": {"sourceDbName": "ChEMBL"}}],
                        },
                        {
                            "drug": {"name": "PLX-4720", "approved": False},
                            "interactionTypes": [
                                {"type": "inhibitor", "directionality": "inhibitory"}
                            ],
                            "interactionClaims": [{"source": {"sourceDbName": "DrugBank"}}],
                        },
                    ],
                }
            ]
        }
    }
}


@respx.mock
async def test_dgidb_get_drug_interactions(http_client):
    respx.post("https://dgidb.org/api/graphql").mock(
        return_value=httpx.Response(200, json=_MOCK_DGIDB_RESPONSE)
    )
    client = DGIdbClient(http_client)
    result = await client.get_drug_interactions("BRAF")

    assert len(result) == 2
    names = {d.drug_name for d in result}
    assert "Vemurafenib" in names
    # Approved drug listed first
    assert result[0].approved is True
    assert result[0].drug_name == "Vemurafenib"
    assert result[0].interaction_type == "inhibitor"
    assert result[0].phase == 4


@respx.mock
async def test_dgidb_returns_empty_on_no_data(http_client):
    respx.post("https://dgidb.org/api/graphql").mock(
        return_value=httpx.Response(200, json={"data": {"genes": {"nodes": []}}})
    )
    client = DGIdbClient(http_client)
    result = await client.get_drug_interactions("FAKEGENE")
    assert result == []


# ---------------------------------------------------------------------------
# ClinicalTrials.gov client tests
# ---------------------------------------------------------------------------


_MOCK_CT_RESPONSE = {
    "studies": [
        {
            "protocolSection": {
                "identificationModule": {
                    "nctId": "NCT01234567",
                    "briefTitle": "Phase 2 Study of Vemurafenib in BRAF V600E Melanoma",
                },
                "statusModule": {"overallStatus": "COMPLETED"},
                "designModule": {"phases": ["PHASE2"]},
                "conditionsModule": {"conditions": ["Melanoma"]},
            }
        },
        {
            "protocolSection": {
                "identificationModule": {
                    "nctId": "NCT09876543",
                    "briefTitle": "Dabrafenib + Trametinib in BRAF Mutant Tumors",
                },
                "statusModule": {"overallStatus": "RECRUITING"},
                "designModule": {"phases": ["PHASE3"]},
                "conditionsModule": {"conditions": ["Non-small Cell Lung Cancer"]},
            }
        },
    ]
}


@respx.mock
async def test_clinical_trials_get_trials(http_client):
    respx.get(url__regex=r"clinicaltrials\.gov/api/v2/studies").mock(
        return_value=httpx.Response(200, json=_MOCK_CT_RESPONSE)
    )
    client = ClinicalTrialsClient(http_client)
    trials, counts = await client.get_trials("BRAF")

    assert len(trials) == 2
    assert trials[0].nct_id == "NCT01234567"
    assert trials[0].phase == "Phase 2"
    assert trials[0].status == "COMPLETED"
    assert "Phase 2" in counts
    assert "Phase 3" in counts
    assert counts["Phase 2"] == 1
    assert counts["Phase 3"] == 1


@respx.mock
async def test_clinical_trials_returns_empty_on_error(http_client):
    respx.get(url__regex=r"clinicaltrials\.gov/api/v2/studies").mock(
        return_value=httpx.Response(500)
    )
    client = ClinicalTrialsClient(http_client)
    trials, counts = await client.get_trials("BRAF")
    assert trials == []
    assert counts == {}


# ---------------------------------------------------------------------------
# Reactome client tests
# ---------------------------------------------------------------------------


_MOCK_REACTOME_TOKEN_RESPONSE = {
    "summary": {
        "token": "MjAyNS0wMS0wMSAxMjo",
        "sampleName": "BRAF",
        "type": "OVERREPRESENTATION",
    },
    "pathways": [
        {
            "stId": "R-HSA-5673001",
            "name": "RAF/MAP kinase cascade",
            "entities": {"pValue": 1.2e-15, "total": 42},
        },
        {
            "stId": "R-HSA-162582",
            "name": "Signal Transduction",
            "entities": {"pValue": 0.0023, "total": 512},
        },
    ],
}


@respx.mock
async def test_reactome_get_pathway_context(http_client):
    # POST to analysis service returns inline pathways — no separate download needed
    respx.post(url__regex=r"reactome\.org/AnalysisService/identifiers").mock(
        return_value=httpx.Response(200, json=_MOCK_REACTOME_TOKEN_RESPONSE)
    )
    client = ReactomeClient(http_client)
    result = await client.get_pathway_context("BRAF")

    assert result is not None
    assert result.gene_symbol == "BRAF"
    assert len(result.pathways) >= 1
    assert result.top_pathway_name is not None
    names = [p.display_name for p in result.pathways]
    assert any("kinase" in n.lower() or "signal" in n.lower() for n in names)


@respx.mock
async def test_reactome_returns_none_on_failure(http_client):
    respx.post(url__regex=r"reactome\.org/AnalysisService/identifiers").mock(
        return_value=httpx.Response(500)
    )
    client = ReactomeClient(http_client)
    result = await client.get_pathway_context("BRAF")
    assert result is None


@respx.mock
async def test_reactome_failure_not_cached(http_client):
    """A transient network error must not poison the cache for subsequent calls."""
    route = respx.post(url__regex=r"reactome\.org/AnalysisService/identifiers")
    # First call fails
    route.mock(return_value=httpx.Response(500))
    client = ReactomeClient(http_client)
    result1 = await client.get_pathway_context("BRAF")
    assert result1 is None

    # Second call succeeds — must not be blocked by a cached None
    route.mock(return_value=httpx.Response(200, json=_MOCK_REACTOME_TOKEN_RESPONSE))
    result2 = await client.get_pathway_context("BRAF")
    assert result2 is not None
    assert result2.gene_symbol == "BRAF"


@respx.mock
async def test_reactome_token_fallback(http_client):
    """When POST inline pathways are empty, fall back to GET /token/{token}/pathways/TOTAL/."""
    _post_no_inline = {
        "summary": {"token": "abc123", "type": "OVERREPRESENTATION"},
        "pathwaysFound": 2,
        "pathways": [],
    }
    _token_pathways = [
        {
            "stId": "R-HSA-5673001",
            "name": "RAF/MAP kinase cascade",
            "entities": {"pValue": 1.2e-15, "total": 42},
        },
    ]
    respx.post(url__regex=r"reactome\.org/AnalysisService/identifiers").mock(
        return_value=httpx.Response(200, json=_post_no_inline)
    )
    respx.get(url__regex=r"reactome\.org/AnalysisService/token/abc123/pathways/TOTAL").mock(
        return_value=httpx.Response(200, json=_token_pathways)
    )
    client = ReactomeClient(http_client)
    result = await client.get_pathway_context("BRAF")

    assert result is not None
    assert result.gene_symbol == "BRAF"
    assert len(result.pathways) == 1
    assert result.pathways[0].reactome_id == "R-HSA-5673001"


@respx.mock
async def test_reactome_token_fallback_no_token_returns_none(http_client):
    """If POST returns empty pathways and no token, return None gracefully."""
    _post_no_token = {
        "summary": {"type": "OVERREPRESENTATION"},
        "pathwaysFound": 0,
        "pathways": [],
    }
    respx.post(url__regex=r"reactome\.org/AnalysisService/identifiers").mock(
        return_value=httpx.Response(200, json=_post_no_token)
    )
    client = ReactomeClient(http_client)
    result = await client.get_pathway_context("UNKNOWNGENE")
    assert result is None


_MOCK_REACTOME_SEARCH_RESPONSE = {
    "results": [
        {
            "entries": [
                {"stId": "R-HSA-5673001", "name": "RAF/MAP kinase cascade"},
            ]
        }
    ]
}

_MOCK_REACTOME_PARTICIPANTS_RESPONSE = [
    {
        "refEntities": [
            {"geneName": ["BRAF"], "className": "ReferenceGeneProduct"},
            {"geneName": ["RAF1"], "className": "ReferenceGeneProduct"},
        ]
    },
    {
        "refEntities": [
            {"geneName": ["MAP2K1", "MEK1"], "className": "ReferenceGeneProduct"},
        ]
    },
]


@respx.mock
async def test_reactome_get_pathway_members_by_name(http_client):
    respx.get(url__regex=r"reactome\.org/ContentService/data/search/query").mock(
        return_value=httpx.Response(200, json=_MOCK_REACTOME_SEARCH_RESPONSE)
    )
    respx.get(url__regex=r"reactome\.org/ContentService/data/participants").mock(
        return_value=httpx.Response(200, json=_MOCK_REACTOME_PARTICIPANTS_RESPONSE)
    )
    client = ReactomeClient(http_client)
    genes = await client.get_pathway_members("RAF/MAP kinase cascade")

    assert "BRAF" in genes
    assert "RAF1" in genes
    # MAP2K1 should appear (first geneName entry wins)
    assert "MAP2K1" in genes


@respx.mock
async def test_reactome_get_pathway_members_by_stid(http_client):
    """Stable ID input skips the search step entirely."""
    respx.get(url__regex=r"reactome\.org/ContentService/data/participants").mock(
        return_value=httpx.Response(200, json=_MOCK_REACTOME_PARTICIPANTS_RESPONSE)
    )
    client = ReactomeClient(http_client)
    genes = await client.get_pathway_members("R-HSA-5673001")

    assert "BRAF" in genes
    assert "RAF1" in genes


@respx.mock
async def test_reactome_get_pathway_members_not_found(http_client):
    respx.get(url__regex=r"reactome\.org/ContentService/data/search/query").mock(
        return_value=httpx.Response(200, json={"results": []})
    )
    client = ReactomeClient(http_client)
    genes = await client.get_pathway_members("nonexistent pathway xyz")
    assert genes == []


# ---------------------------------------------------------------------------
# EFOResolver tests
# ---------------------------------------------------------------------------

_MOCK_OLS4_OBESITY = {
    "response": {
        "docs": [
            {
                "iri": "http://www.ebi.ac.uk/efo/EFO_0001073",
                "label": "obesity",
                "synonym": ["overweight", "adiposity", "fat body mass"],
            },
            {
                "iri": "http://www.ebi.ac.uk/efo/EFO_0004340",
                "label": "body mass index",
                "synonym": ["BMI"],
            },
        ]
    }
}


@respx.mock
async def test_efo_resolver_returns_terms(http_client):
    respx.get(url__regex=r"ols4/api/search").mock(
        return_value=httpx.Response(200, json=_MOCK_OLS4_OBESITY)
    )
    resolver = EFOResolver(http_client, cache_path=None)
    terms = await resolver.resolve("obesity")

    assert len(terms) == 2
    assert terms[0].uri == "http://www.ebi.ac.uk/efo/EFO_0001073"
    assert terms[0].label == "obesity"
    assert "adiposity" in terms[0].synonyms
    assert terms[1].label == "body mass index"


@respx.mock
async def test_efo_resolver_falls_back_on_ols_failure(http_client):
    respx.get(url__regex=r"ols4/api/search").mock(return_value=httpx.Response(503))
    resolver = EFOResolver(http_client, cache_path=None)
    terms = await resolver.resolve("obesity")
    # Must return [] gracefully — callers fall back to synonym dict
    assert terms == []


@respx.mock
async def test_efo_resolver_session_cache_prevents_duplicate_calls(http_client):
    call_count = 0

    def side_effect(request):
        nonlocal call_count
        call_count += 1
        return httpx.Response(200, json=_MOCK_OLS4_OBESITY)

    respx.get(url__regex=r"ols4/api/search").mock(side_effect=side_effect)
    resolver = EFOResolver(http_client, cache_path=None)

    await resolver.resolve("obesity")
    await resolver.resolve("obesity")  # second call — must hit session cache

    assert call_count == 1


# ---------------------------------------------------------------------------
# filter_by_trait tests
# ---------------------------------------------------------------------------


def _make_hit(trait: str, efo_uri: str | None = None) -> GwasHit:
    return GwasHit(
        study_accession="GCST000001",
        trait=trait,
        mapped_gene="FTO",
        risk_allele="rs1234-A",
        p_value=1e-10,
        efo_uri=efo_uri,
    )


def test_filter_by_trait_efo_uri_exact_match():
    """EFO URI hit matches precisely — string content of trait label is irrelevant."""
    hits = [
        _make_hit("body mass index", efo_uri="http://www.ebi.ac.uk/efo/EFO_0001073"),
        _make_hit("type 2 diabetes", efo_uri="http://www.ebi.ac.uk/efo/EFO_0000400"),
    ]
    efo_terms = [EFOTerm(uri="http://www.ebi.ac.uk/efo/EFO_0001073", label="obesity", synonyms=[])]
    result = filter_by_trait(hits, "fat", efo_terms=efo_terms)

    assert len(result) == 1
    assert result[0].trait == "body mass index"


def test_filter_by_trait_efo_synonym_string_match():
    """EFO synonym expansion catches label variants not in the hardcoded dict."""
    hits = [
        _make_hit("fat body mass"),  # no efo_uri — SNP path
        _make_hit("coronary artery disease"),
    ]
    efo_terms = [
        EFOTerm(
            uri="http://www.ebi.ac.uk/efo/EFO_0001073",
            label="obesity",
            synonyms=["fat body mass", "adiposity"],
        )
    ]
    result = filter_by_trait(hits, "fat", efo_terms=efo_terms)

    assert len(result) == 1
    assert result[0].trait == "fat body mass"


def test_filter_by_trait_fallback_to_synonym_dict_when_no_efo():
    """When efo_terms is None, falls back to hardcoded TRAIT_SYNONYMS dict."""
    hits = [
        _make_hit("low-density lipoprotein cholesterol"),
        _make_hit("type 2 diabetes"),
    ]
    result = filter_by_trait(hits, "hypercholesterolemia", efo_terms=None)

    assert len(result) == 1
    assert "cholesterol" in result[0].trait


def test_filter_by_trait_unknown_indication_direct_substring():
    """Indication not in dict and no EFO terms — falls back to direct substring match."""
    hits = [
        _make_hit("autoimmune thyroid disease"),
        _make_hit("melanoma"),
    ]
    result = filter_by_trait(hits, "thyroid disease", efo_terms=None)

    assert len(result) == 1
    assert result[0].trait == "autoimmune thyroid disease"


# ---------------------------------------------------------------------------
# SAbDab client tests
# ---------------------------------------------------------------------------

_MOCK_SABDAB_TSV = (
    "pdb\tHchain\tLchain\tmodel\tantigen_chain\tantigen_type\tantigen_het_name\t"
    "antigen_name\tshort_header\tdate\tcompound\torganism\theavy_species\tlight_species\t"
    "antigen_species\tauthors\tresolution\tmethod\tr_free\tr_factor\tscfv\tengineered\t"
    "heavy_subclass\tlight_subclass\tlight_ctype\taffinity\tdelta_g\taffinity_method\t"
    "temperature\tpmid\n"
    "1abc\tH\tL\t0\tA\tprotein\tNA\tepidermal growth factor receptor\tIMMUNE SYSTEM\t"
    "01/01/20\tAnti-EGFR antibody complex\tHomo sapiens\thomo sapiens\thomo sapiens\t"
    "homo sapiens\tSmith, J.\t2.5\tX-RAY DIFFRACTION\t0.21\t0.18\tFalse\tTrue\t"
    "IGHV3\tIGKV1\tKappa\t5.2\tNone\tITC\t25\t12345678\n"
    "2def\tH\tNA\t0\tA\tprotein\tNA\tepidermal growth factor receptor\tIMMUNE SYSTEM\t"
    "06/15/21\tEGFR nanobody VHH complex\tLama glama\tlama glama\t\t"
    "homo sapiens\tJones, A.\t3.1\tELECTRON MICROSCOPY\tNA\tNA\tFalse\tTrue\t"
    "IGHV1\tNA\tNA\tNone\tNone\tNone\tNone\tNone\n"
    "3ghi\tH\tL\t0\tB\tprotein\tNA\ttumor necrosis factor\tIMMUNE SYSTEM\t"
    "03/10/22\tAnti-TNF Fab fragment\tHomo sapiens\thomo sapiens\thomo sapiens\t"
    "homo sapiens\tBrown, K.\t1.9\tX-RAY DIFFRACTION\t0.19\t0.16\tFalse\tTrue\t"
    "IGHV1\tIGLV2\tLambda\tNone\tNone\tNone\tNone\t99999999\n"
)


@respx.mock
@pytest.mark.asyncio
async def test_sabdab_get_antibody_structures_happy_path(http_client, tmp_path, monkeypatch):
    """Should parse TSV, filter by antigen name, and return AntibodyStructures."""
    monkeypatch.setattr(
        "genesis_bio_mcp.clients.sabdab.settings.sabdab_cache_path",
        tmp_path / "sabdab_cache.tsv",
    )
    monkeypatch.setattr(
        "genesis_bio_mcp.clients.sabdab.settings.sabdab_cache_ttl_secs",
        604800,
    )

    respx.get(url__regex=r"sabdab-sabpred.*summary").mock(
        return_value=httpx.Response(200, content=_MOCK_SABDAB_TSV.encode())
    )

    client = SAbDabClient(http_client)
    result = await client.get_antibody_structures("egfr")

    assert result is not None
    assert result.total_structures == 2
    assert result.nanobody_count == 1
    assert result.fab_count == 1
    # Best resolution (2.5 Å X-ray) should be first
    assert result.structures[0].pdb == "1ABC"
    assert result.structures[0].is_nanobody is False
    assert result.structures[0].resolution_ang == pytest.approx(2.5)
    assert result.structures[0].is_engineered is True
    assert result.structures[0].affinity_nM == pytest.approx(5.2)
    # Second is the VHH
    assert result.structures[1].pdb == "2DEF"
    assert result.structures[1].is_nanobody is True


@respx.mock
@pytest.mark.asyncio
async def test_sabdab_empty_results(http_client, tmp_path, monkeypatch):
    """Returns AntibodyStructures with total_structures=0 when no match."""
    monkeypatch.setattr(
        "genesis_bio_mcp.clients.sabdab.settings.sabdab_cache_path",
        tmp_path / "sabdab_cache.tsv",
    )
    monkeypatch.setattr(
        "genesis_bio_mcp.clients.sabdab.settings.sabdab_cache_ttl_secs",
        604800,
    )

    respx.get(url__regex=r"sabdab-sabpred.*summary").mock(
        return_value=httpx.Response(200, content=_MOCK_SABDAB_TSV.encode())
    )

    client = SAbDabClient(http_client)
    result = await client.get_antibody_structures("nonexistent_target_xyz")

    assert result is not None
    assert result.total_structures == 0
    assert result.structures == []


@respx.mock
@pytest.mark.asyncio
async def test_sabdab_download_failure_returns_none(http_client, tmp_path, monkeypatch):
    """Returns None when download fails and no disk cache exists."""
    monkeypatch.setattr(
        "genesis_bio_mcp.clients.sabdab.settings.sabdab_cache_path",
        tmp_path / "sabdab_cache.tsv",
    )
    monkeypatch.setattr(
        "genesis_bio_mcp.clients.sabdab.settings.sabdab_cache_ttl_secs",
        604800,
    )

    respx.get(url__regex=r"sabdab-sabpred.*summary").mock(return_value=httpx.Response(500))

    client = SAbDabClient(http_client)
    result = await client.get_antibody_structures("egfr")

    assert result is None


@pytest.mark.asyncio
async def test_sabdab_uses_disk_cache(http_client, tmp_path, monkeypatch):
    """If a fresh disk cache exists, no HTTP request is made."""
    cache_path = tmp_path / "sabdab_cache.tsv"
    cache_path.write_bytes(_MOCK_SABDAB_TSV.encode())

    monkeypatch.setattr(
        "genesis_bio_mcp.clients.sabdab.settings.sabdab_cache_path",
        cache_path,
    )
    monkeypatch.setattr(
        "genesis_bio_mcp.clients.sabdab.settings.sabdab_cache_ttl_secs",
        604800,
    )

    client = SAbDabClient(http_client)
    result = await client.get_antibody_structures("tnf")

    assert result is not None
    assert result.total_structures == 1
    assert result.structures[0].pdb == "3GHI"


def test_sabdab_to_markdown():
    """AntibodyStructures.to_markdown returns a non-empty string with expected content."""
    from genesis_bio_mcp.models import AntibodyStructure, AntibodyStructures

    structures = [
        AntibodyStructure(
            pdb="1ABC",
            is_nanobody=False,
            antigen_name="epidermal growth factor receptor",
            resolution_ang=2.5,
            method="X-RAY DIFFRACTION",
            heavy_species="homo sapiens",
            light_species="homo sapiens",
            heavy_subclass="IGHV3",
            light_subclass="IGKV1",
            is_engineered=True,
            is_scfv=False,
            affinity_nM=5.2,
            compound="Anti-EGFR Fab",
            date_added="01/01/20",
            pmid="12345678",
        ),
        AntibodyStructure(
            pdb="2DEF",
            is_nanobody=True,
            antigen_name="epidermal growth factor receptor",
            resolution_ang=3.1,
            method="ELECTRON MICROSCOPY",
            heavy_species="lama glama",
            light_species=None,
            heavy_subclass="IGHV1",
            light_subclass=None,
            is_engineered=True,
            is_scfv=False,
            affinity_nM=None,
            compound="EGFR VHH nanobody",
            date_added="06/15/21",
            pmid=None,
        ),
    ]
    result = AntibodyStructures(
        query="EGFR",
        total_structures=2,
        nanobody_count=1,
        fab_count=1,
        structures=structures,
    )
    md = result.to_markdown()
    assert "EGFR" in md
    assert "1ABC" in md
    assert "2DEF" in md
    assert "VHH" in md
    assert "2.50 Å" in md


# Minimal RCSB FASTA for 1ABC: chain H (heavy) and chain L (light)
_MOCK_RCSB_FASTA = (
    ">1ABC_1|Chains H|Anti-EGFR heavy chain|Homo sapiens (9606)\n"
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVK"
    "GRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS\n"
    ">1ABC_2|Chains L|Anti-EGFR light chain|Homo sapiens (9606)\n"
    "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLSWYQQKPGKAPKLLIYAASSLQSGVPSRFSG"
    "SRSGTDFTLTISSLQPEDFATYYCQQSYSTPPTFGQGTKVEIK\n"
)

# Minimal AbNum output (Chothia) for heavy chain — only CDR-relevant positions
_MOCK_ABNUM_VH = "\n".join(
    [
        # CDR-H1 (Chothia H26-H32): GFNIKDT
        "H26 G",
        "H27 F",
        "H28 N",
        "H29 I",
        "H30 K",
        "H31 D",
        "H32 T",
        # CDR-H2 (Chothia H52-H58): IYPTNG + insertions
        "H52 I",
        "H52A Y",
        "H53 P",
        "H54 T",
        "H55 N",
        "H56 G",
        "H57 Y",
        "H58 T",
        # CDR-H3 (Chothia H95-H102): RWGGDGFY + insertion
        "H95 R",
        "H96 W",
        "H97 G",
        "H98 G",
        "H99 D",
        "H100 G",
        "H100A F",
        "H101 Y",
        "H102 A",
    ]
)

# Minimal AbNum output for light chain (Chothia): CDR-L1/2/3
_MOCK_ABNUM_VL = "\n".join(
    [
        # CDR-L1 (Chothia L24-L34): RASQSISSYLS
        "L24 R",
        "L25 A",
        "L26 S",
        "L27 Q",
        "L28 S",
        "L29 I",
        "L30 S",
        "L31 S",
        "L32 Y",
        "L33 L",
        "L34 S",
        # CDR-L2 (Chothia L50-L56): AASSLQS
        "L50 A",
        "L51 A",
        "L52 S",
        "L53 S",
        "L54 L",
        "L55 Q",
        "L56 S",
        # CDR-L3 (Chothia L89-L97): QQSYSTPPT
        "L89 Q",
        "L90 Q",
        "L91 S",
        "L92 Y",
        "L93 S",
        "L94 T",
        "L95 P",
        "L96 P",
        "L97 T",
    ]
)


@respx.mock
@pytest.mark.asyncio
async def test_sabdab_cdr_annotation_happy_path(http_client, tmp_path, monkeypatch):
    """Top structures should have VH and VL CDR sequences populated via AbNum."""
    monkeypatch.setattr(
        "genesis_bio_mcp.clients.sabdab.settings.sabdab_cache_path",
        tmp_path / "sabdab_cache.tsv",
    )
    monkeypatch.setattr(
        "genesis_bio_mcp.clients.sabdab.settings.sabdab_cache_ttl_secs",
        604800,
    )

    respx.get(url__regex=r"sabdab-sabpred.*summary").mock(
        return_value=httpx.Response(200, content=_MOCK_SABDAB_TSV.encode())
    )
    respx.get(url__regex=r"rcsb\.org/fasta").mock(
        return_value=httpx.Response(200, text=_MOCK_RCSB_FASTA)
    )

    # Mock AbNum: return VH output for H chain, VL for L chain
    def _abnum_side_effect(request, **kwargs):
        aaseq = dict(request.url.params).get("aaseq", "")
        # Identify chain type by sequence content (crude but deterministic in tests)
        return httpx.Response(200, text=_MOCK_ABNUM_VH if "EVQL" in aaseq else _MOCK_ABNUM_VL)

    respx.get(url__regex=r"bioinf\.org\.uk.*abnum").mock(side_effect=_abnum_side_effect)

    client = SAbDabClient(http_client)
    result = await client.get_antibody_structures("egfr")

    assert result is not None
    top = result.structures[0]
    assert top.pdb == "1ABC"
    assert top.vh_cdr1 == "GFNIKDT"
    assert top.vh_cdr2 == "IYP TNGYT".replace(" ", "")  # H52-H58 with H52A
    assert top.vh_cdr3 == "RWGGDGFYA"  # H95-H102 with H100A
    assert top.vl_cdr1 == "RASQSISSYLS"
    assert top.vl_cdr2 == "AASSLQS"
    assert top.vl_cdr3 == "QQSYSTPPT"


@respx.mock
@pytest.mark.asyncio
async def test_sabdab_cdr_annotation_fasta_failure(http_client, tmp_path, monkeypatch):
    """When RCSB FASTA fails, structures are returned with CDR fields as None."""
    monkeypatch.setattr(
        "genesis_bio_mcp.clients.sabdab.settings.sabdab_cache_path",
        tmp_path / "sabdab_cache.tsv",
    )
    monkeypatch.setattr(
        "genesis_bio_mcp.clients.sabdab.settings.sabdab_cache_ttl_secs",
        604800,
    )

    respx.get(url__regex=r"sabdab-sabpred.*summary").mock(
        return_value=httpx.Response(200, content=_MOCK_SABDAB_TSV.encode())
    )
    respx.get(url__regex=r"rcsb\.org/fasta").mock(return_value=httpx.Response(404))

    client = SAbDabClient(http_client)
    result = await client.get_antibody_structures("egfr")

    assert result is not None
    assert result.total_structures == 2
    top = result.structures[0]
    assert top.pdb == "1ABC"
    assert top.vh_cdr1 is None
    assert top.vh_cdr2 is None
    assert top.vh_cdr3 is None


@respx.mock
@pytest.mark.asyncio
async def test_sabdab_cdr_annotation_abnum_failure(http_client, tmp_path, monkeypatch):
    """When AbNum returns empty output, CDR fields are None and no exception propagates."""
    monkeypatch.setattr(
        "genesis_bio_mcp.clients.sabdab.settings.sabdab_cache_path",
        tmp_path / "sabdab_cache.tsv",
    )
    monkeypatch.setattr(
        "genesis_bio_mcp.clients.sabdab.settings.sabdab_cache_ttl_secs",
        604800,
    )

    respx.get(url__regex=r"sabdab-sabpred.*summary").mock(
        return_value=httpx.Response(200, content=_MOCK_SABDAB_TSV.encode())
    )
    respx.get(url__regex=r"rcsb\.org/fasta").mock(
        return_value=httpx.Response(200, text=_MOCK_RCSB_FASTA)
    )
    respx.get(url__regex=r"bioinf\.org\.uk.*abnum").mock(
        return_value=httpx.Response(200, text="")  # empty — no numbering assigned
    )

    client = SAbDabClient(http_client)
    result = await client.get_antibody_structures("egfr")

    assert result is not None
    top = result.structures[0]
    assert top.vh_cdr1 is None
    assert top.vh_cdr2 is None
    assert top.vh_cdr3 is None
    assert top.vl_cdr1 is None


# ---------------------------------------------------------------------------
# gnomAD client tests
# ---------------------------------------------------------------------------

_MOCK_GNOMAD_RESPONSE = {
    "data": {
        "gene": {
            "gene_id": "ENSG00000157764",
            "name": "B-Raf proto-oncogene, serine/threonine kinase",
            "canonical_transcript_id": "ENST00000644969",
            "gnomad_constraint": {
                "pLI": 0.9999,
                "lof_z": 7.35,
                "mis_z": 5.52,
                "oe_lof": 0.153,
                "oe_lof_lower": 0.104,
                "oe_lof_upper": 0.232,
                "oe_mis": 0.577,
                "exp_lof": 104.6,
                "exp_mis": 937.3,
                "obs_lof": 16,
                "obs_mis": 541,
            },
        }
    }
}

_MOCK_GNOMAD_NO_CONSTRAINT = {
    "data": {
        "gene": {
            "gene_id": "ENSG00000123456",
            "name": "Hypothetical gene",
            "canonical_transcript_id": None,
            "gnomad_constraint": None,
        }
    }
}


@respx.mock
@pytest.mark.asyncio
async def test_gnomad_get_constraint_happy_path(http_client):
    """Should parse gnomAD constraint response and return GnomADConstraint."""
    respx.post(url__regex=r"gnomad\.broadinstitute\.org/api").mock(
        return_value=httpx.Response(200, json=_MOCK_GNOMAD_RESPONSE)
    )
    client = GnomADClient(http_client)
    result = await client.get_constraint("BRAF")

    assert result is not None
    assert result.constraint_available is True
    assert result.gene_symbol == "BRAF"
    assert result.ensembl_id == "ENSG00000157764"
    assert result.pLI == pytest.approx(0.9999)
    assert result.oe_lof_upper == pytest.approx(0.232)
    assert result.obs_lof == 16


@respx.mock
@pytest.mark.asyncio
async def test_gnomad_no_constraint_data(http_client):
    """Returns GnomADConstraint with constraint_available=False when gnomad_constraint is null."""
    respx.post(url__regex=r"gnomad\.broadinstitute\.org/api").mock(
        return_value=httpx.Response(200, json=_MOCK_GNOMAD_NO_CONSTRAINT)
    )
    client = GnomADClient(http_client)
    result = await client.get_constraint("FAKEGENE")

    assert result is not None
    assert result.constraint_available is False
    assert result.pLI is None


@respx.mock
@pytest.mark.asyncio
async def test_gnomad_network_error_returns_none(http_client):
    """Returns None on network error."""
    respx.post(url__regex=r"gnomad\.broadinstitute\.org/api").mock(
        side_effect=httpx.ConnectError("timeout")
    )
    client = GnomADClient(http_client)
    result = await client.get_constraint("BRAF")

    assert result is None


def test_gnomad_to_markdown_constrained():
    """GnomADConstraint.to_markdown includes pLI, LOEUF, and engineering note."""
    from genesis_bio_mcp.models import GnomADConstraint

    c = GnomADConstraint(
        gene_symbol="BRAF",
        ensembl_id="ENSG00000157764",
        gene_name="B-Raf proto-oncogene",
        constraint_available=True,
        pLI=0.9999,
        lof_z=7.35,
        mis_z=5.52,
        oe_lof=0.153,
        oe_lof_lower=0.104,
        oe_lof_upper=0.232,
        oe_mis=0.577,
        exp_lof=104.6,
        exp_mis=937.3,
        obs_lof=16,
        obs_mis=541,
    )
    md = c.to_markdown()
    assert "BRAF" in md
    assert "pLI" in md
    assert "LOEUF" in md
    assert "Highly constrained" in md
    assert "Engineering note" in md


# ---------------------------------------------------------------------------
# InterPro client tests
# ---------------------------------------------------------------------------

_MOCK_INTERPRO_RESPONSE = {
    "count": 2,
    "next": None,
    "previous": None,
    "results": [
        {
            "metadata": {
                "accession": "IPR000719",
                "name": "Protein kinase domain",
                "source_database": "interpro",
                "type": "domain",
                "member_databases": {
                    "pfam": {"PF00069": "Protein kinase domain"},
                    "smart": {"SM00220": "S_TKc"},
                },
                "go_terms": [
                    {
                        "identifier": "GO:0004672",
                        "name": "protein kinase activity",
                        "category": {"code": "F", "name": "molecular_function"},
                    },
                ],
            },
            "proteins": [
                {
                    "accession": "p15056",
                    "protein_length": 766,
                    "source_database": "reviewed",
                    "organism": "9606",
                    "entry_protein_locations": [
                        {
                            "fragments": [{"start": 457, "end": 717, "dc-status": "CONTINUOUS"}],
                            "representative": True,
                            "model": "IPR000719",
                            "score": 0,
                        }
                    ],
                }
            ],
        },
        {
            "metadata": {
                "accession": "IPR003116",
                "name": "Raf-like Ras-binding",
                "source_database": "interpro",
                "type": "domain",
                "member_databases": {"pfam": {"PF02196": "RBD"}},
                "go_terms": [],
            },
            "proteins": [
                {
                    "accession": "p15056",
                    "protein_length": 766,
                    "source_database": "reviewed",
                    "organism": "9606",
                    "entry_protein_locations": [
                        {
                            "fragments": [{"start": 155, "end": 227, "dc-status": "CONTINUOUS"}],
                            "representative": True,
                            "model": "IPR003116",
                            "score": 0,
                        }
                    ],
                }
            ],
        },
    ],
}


@respx.mock
@pytest.mark.asyncio
async def test_interpro_get_domains_happy_path(http_client):
    """Should parse InterPro response and return DomainAnnotations sorted by position."""
    respx.get(url__regex=r"ebi\.ac\.uk/interpro/api/entry/InterPro/protein").mock(
        return_value=httpx.Response(200, json=_MOCK_INTERPRO_RESPONSE)
    )
    client = InterProClient(http_client)
    result = await client.get_domains("BRAF", "P15056")

    assert result is not None
    assert result.total_entries == 2
    assert len(result.domains) == 2
    # Should be sorted by start position: RBD (155) before kinase (457)
    assert result.domains[0].interpro_accession == "IPR003116"
    assert result.domains[1].interpro_accession == "IPR000719"
    assert result.domains[1].positions == [(457, 717)]
    assert "PF00069" in result.domains[1].member_databases.get("pfam", [])


@respx.mock
@pytest.mark.asyncio
async def test_interpro_returns_empty_on_404(http_client):
    """Returns DomainAnnotations with 0 entries on 404."""
    respx.get(url__regex=r"ebi\.ac\.uk/interpro/api/entry/InterPro/protein").mock(
        return_value=httpx.Response(404)
    )
    client = InterProClient(http_client)
    result = await client.get_domains("FAKEGENE", "X00000")

    assert result is not None
    assert result.total_entries == 0
    assert result.domains == []


@respx.mock
@pytest.mark.asyncio
async def test_interpro_returns_none_on_error(http_client):
    """Returns None on network error."""
    respx.get(url__regex=r"ebi\.ac\.uk/interpro/api/entry/InterPro/protein").mock(
        side_effect=httpx.ConnectError("timeout")
    )
    client = InterProClient(http_client)
    result = await client.get_domains("BRAF", "P15056")

    assert result is None


def test_interpro_to_markdown():
    """DomainAnnotations.to_markdown includes domain names and positions."""
    from genesis_bio_mcp.models import DomainAnnotation, DomainAnnotations

    result = DomainAnnotations(
        gene_symbol="BRAF",
        uniprot_accession="P15056",
        total_entries=1,
        domains=[
            DomainAnnotation(
                interpro_accession="IPR000719",
                name="Protein kinase domain",
                entry_type="domain",
                positions=[(457, 717)],
                member_databases={"pfam": ["PF00069"]},
                go_terms=["GO:0004672 protein kinase activity"],
            )
        ],
    )
    md = result.to_markdown()
    assert "BRAF" in md
    assert "IPR000719" in md
    assert "457" in md
    assert "717" in md
    assert "PF00069" in md


# ---------------------------------------------------------------------------
# IEDB client tests
# ---------------------------------------------------------------------------

_MOCK_IEDB_RESPONSE = [
    {
        "structure_description": "K443, L444, I448, M480",
        "linear_sequence": None,
        "qualitative_measure": "Positive",
        "antibody_isotype": "IgG1",
        "pubmed_id": "15837620",
        "pdb_id": "1yy9",
        "curated_source_antigen": {
            "accession": "P00533.2",
            "name": "Epidermal growth factor receptor",
            "starting_position": None,
            "ending_position": None,
        },
    },
    {
        "structure_description": "EEEEEEIVYK",
        "linear_sequence": "EEEEEEIVYK",
        "qualitative_measure": "Positive",
        "antibody_isotype": "IgG",
        "pubmed_id": "12345678",
        "pdb_id": None,
        "curated_source_antigen": {
            "accession": "P00533.2",
            "name": "Epidermal growth factor receptor",
            "starting_position": 100,
            "ending_position": 109,
        },
    },
]


@respx.mock
@pytest.mark.asyncio
async def test_iedb_get_epitopes_happy_path(http_client):
    """Should parse IEDB response and return EpitopeResults."""
    respx.get(url__regex=r"query-api\.iedb\.org/api/v1/bcell_search").mock(
        return_value=httpx.Response(200, json=_MOCK_IEDB_RESPONSE)
    )
    client = IEDBClient(http_client)
    result = await client.get_epitopes("epidermal growth factor receptor")

    assert result is not None
    assert result.total_assays == 2
    assert result.unique_epitopes == 2
    assert result.with_structure == 1
    assert len(result.epitopes) == 2


@respx.mock
@pytest.mark.asyncio
async def test_iedb_empty_results(http_client):
    """Returns EpitopeResults with 0 assays on empty list response."""
    respx.get(url__regex=r"query-api\.iedb\.org/api/v1/bcell_search").mock(
        return_value=httpx.Response(200, json=[])
    )
    client = IEDBClient(http_client)
    result = await client.get_epitopes("unknown_antigen")

    assert result is not None
    assert result.total_assays == 0
    assert result.epitopes == []


@respx.mock
@pytest.mark.asyncio
async def test_iedb_network_error_returns_none(http_client):
    """Returns None on network error."""
    respx.get(url__regex=r"query-api\.iedb\.org/api/v1/bcell_search").mock(
        side_effect=httpx.ConnectError("timeout")
    )
    client = IEDBClient(http_client)
    result = await client.get_epitopes("egfr")

    assert result is None


def test_iedb_to_markdown():
    """EpitopeResults.to_markdown includes assay count and epitope sequences."""
    from genesis_bio_mcp.models import EpitopeRecord, EpitopeResults

    result = EpitopeResults(
        antigen_query="epidermal growth factor receptor",
        total_assays=2,
        unique_epitopes=2,
        with_structure=1,
        epitopes=[
            EpitopeRecord(
                sequence="K443, L444, I448",
                isotype="IgG1",
                pmid="15837620",
                pdb_id="1yy9",
                antigen_name="Epidermal growth factor receptor",
                antigen_accession="P00533.2",
                start_position=None,
                end_position=None,
            ),
        ],
    )
    md = result.to_markdown()
    assert "epidermal growth factor receptor" in md
    assert "2 positive assays" in md
    assert "1yy9" in md
    assert "IgG1" in md


# ---------------------------------------------------------------------------
# MaveDB client tests
# ---------------------------------------------------------------------------

_MOCK_MAVEDB_SEARCH = {
    "scoreSets": [
        {
            "urn": "urn:mavedb:00000097-a-1",
            "title": "BRCA1 RING domain saturation genome editing",
            "shortDescription": "SGE fitness scores for BRCA1 RING domain variants",
            "numVariants": 3893,
            "publishedDate": "2018-09-12",
            "targetGenes": [
                {
                    "name": "BRCA1",
                    "uniprotIdFromMappedMetadata": "P38398",
                }
            ],
            "primaryPublicationIdentifiers": [
                {"dbName": "PubMed", "identifier": "30209399"},
            ],
        },
        {
            "urn": "urn:mavedb:00000050-a-1",
            "title": "BRCA1 BRCT domain functional scores",
            "shortDescription": "Deep mutational scan of BRCA1 BRCT domain",
            "numVariants": 1056,
            "publishedDate": "2018-01-25",
            "targetGenes": [
                {
                    "name": "BRCA1",
                    "uniprotIdFromMappedMetadata": "P38398",
                }
            ],
            "primaryPublicationIdentifiers": [
                {"dbName": "PubMed", "identifier": "29785012"},
                {"dbName": "doi", "identifier": "10.1016/j.cell.2018.03.072"},
            ],
        },
    ]
}


@respx.mock
@pytest.mark.asyncio
async def test_mavedb_happy_path(http_client):
    """Returns DMSResults with sorted score sets on success."""
    respx.post(url__regex=r"api\.mavedb\.org/api/v1/score-sets/search").mock(
        return_value=httpx.Response(200, json=_MOCK_MAVEDB_SEARCH)
    )
    client = MaveDBClient(http_client)
    result = await client.get_dms_scores("BRCA1")

    assert result is not None
    assert result.gene_symbol == "BRCA1"
    assert result.total_score_sets == 2
    assert result.total_variants == 3893 + 1056
    # Sorted by variant count descending
    assert result.score_sets[0].num_variants == 3893
    assert result.score_sets[0].urn == "urn:mavedb:00000097-a-1"
    assert result.score_sets[0].uniprot_accession == "P38398"
    assert result.score_sets[0].pmid == "30209399"
    assert result.score_sets[1].doi == "10.1016/j.cell.2018.03.072"


@respx.mock
@pytest.mark.asyncio
async def test_mavedb_empty_results(http_client):
    """Returns DMSResults with zero score sets when no data found."""
    respx.post(url__regex=r"api\.mavedb\.org/api/v1/score-sets/search").mock(
        return_value=httpx.Response(200, json={"scoreSets": []})
    )
    client = MaveDBClient(http_client)
    result = await client.get_dms_scores("UNKNOWNGENE")

    assert result is not None
    assert result.total_score_sets == 0
    assert result.score_sets == []


@respx.mock
@pytest.mark.asyncio
async def test_mavedb_network_error_returns_none(http_client):
    """Returns None on network error."""
    respx.post(url__regex=r"api\.mavedb\.org/api/v1/score-sets/search").mock(
        side_effect=httpx.ConnectError("timeout")
    )
    client = MaveDBClient(http_client)
    result = await client.get_dms_scores("BRCA1")

    assert result is None


def test_mavedb_to_markdown():
    """DMSResults.to_markdown includes gene name, counts, and URNs."""
    from genesis_bio_mcp.models import DMSResults, DMSScoreSet

    result = DMSResults(
        gene_symbol="BRCA1",
        total_score_sets=1,
        total_variants=3893,
        score_sets=[
            DMSScoreSet(
                urn="urn:mavedb:00000097-a-1",
                title="BRCA1 RING domain saturation genome editing",
                short_description=None,
                num_variants=3893,
                target_gene="BRCA1",
                uniprot_accession="P38398",
                published_date="2018-09-12",
                pmid="30209399",
                doi=None,
            )
        ],
    )
    md = result.to_markdown()
    assert "BRCA1" in md
    assert "1 score set" in md
    assert "3,893" in md
    assert "urn:mavedb:00000097-a-1" in md
    assert "30209399" in md


# ---------------------------------------------------------------------------
# MyVariant.info client tests
# ---------------------------------------------------------------------------

MOCK_MYVARIANT_TP53_R175H = {
    "_id": "chr17:g.7675088C>T",
    "clinvar": {
        "rsid": "rs28934578",
        "variant_id": 12362,
        "hgvs": {
            "coding": ["NM_000546.6:c.524G>A"],
            "genomic": ["NC_000017.11:g.7675088C>T"],
            "protein": ["NP_000537.3:p.Arg175His", "P04637:p.Arg175His"],
        },
        "rcv": [
            {
                "accession": "RCV000013173",
                "clinical_significance": "Pathogenic",
                "review_status": "criteria provided, multiple submitters, no conflicts",
                "origin": "germline",
                "last_evaluated": "2024-06-12",
                "conditions": {"name": "Li-Fraumeni syndrome 1"},
            },
            {
                "accession": "RCV000204931",
                "clinical_significance": "Pathogenic",
                "review_status": "reviewed by expert panel",
                "origin": "germline",
                "last_evaluated": "2024-09-06",
                "conditions": {"name": "Li-Fraumeni syndrome"},
            },
        ],
    },
    "dbnsfp": {
        "alphamissense": {
            "pred": ["P", "P", "P", "P"],
            "rankscore": 0.924,
            "score": [0.978, 0.984, 0.968, 0.985],
        },
        "revel": {"rankscore": 0.980, "score": [0.922, 0.922]},
        "cadd": {"phred": 25.9, "raw_score": 4.598},
        "sift": {"score": [0.08, 0.07]},
        "polyphen2": {"score": [0.704]},
    },
    "gnomad_exome": {"af": {"af": 3.98e-6, "af_nfe": 8.80e-6, "af_afr": 0.0}},
}


@respx.mock
async def test_myvariant_parses_clinvar_and_alphamissense(http_client):
    respx.get("https://myvariant.info/v1/variant/chr17:g.7675088C>T").mock(
        return_value=httpx.Response(200, json=MOCK_MYVARIANT_TP53_R175H)
    )
    client = MyVariantClient(http_client)
    result = await client.get_annotation("chr17:g.7675088C>T")
    assert result is not None
    assert result.clinvar is not None
    assert result.clinvar.rsid == "rs28934578"
    assert len(result.clinvar.assertions) == 2
    assert result.clinvar.significance_summary == "Pathogenic"
    assert "Li-Fraumeni syndrome 1" in result.clinvar.assertions[0].conditions
    assert result.in_silico is not None
    assert result.in_silico.alphamissense_score == pytest.approx(0.97875, abs=0.001)
    assert result.in_silico.alphamissense_class == "likely_pathogenic"
    assert result.in_silico.revel_score == pytest.approx(0.922, abs=0.001)
    assert result.in_silico.cadd_phred == 25.9
    assert result.gnomad is not None
    assert result.gnomad.overall_af == pytest.approx(3.98e-6)
    assert "af_nfe" in result.gnomad.by_population


@respx.mock
async def test_myvariant_returns_none_on_404(http_client):
    respx.get("https://myvariant.info/v1/variant/chr99:g.1A>T").mock(
        return_value=httpx.Response(404, text="Not found")
    )
    client = MyVariantClient(http_client)
    result = await client.get_annotation("chr99:g.1A>T")
    assert result is None


@respx.mock
async def test_myvariant_returns_none_on_network_error(http_client):
    respx.get("https://myvariant.info/v1/variant/chr17:g.7675088C>T").mock(
        side_effect=httpx.ConnectError("boom")
    )
    client = MyVariantClient(http_client)
    result = await client.get_annotation("chr17:g.7675088C>T")
    assert result is None


@respx.mock
async def test_myvariant_handles_missing_fields(http_client):
    respx.get("https://myvariant.info/v1/variant/chr1:g.1A>T").mock(
        return_value=httpx.Response(200, json={"_id": "chr1:g.1A>T"})
    )
    client = MyVariantClient(http_client)
    result = await client.get_annotation("chr1:g.1A>T")
    # Empty payload → all fields None, no crash
    assert result is not None
    assert result.clinvar is None
    assert result.in_silico is None
    assert result.gnomad is None


# ---------------------------------------------------------------------------
# gnomAD find_variant_id_by_protein_change tests
# ---------------------------------------------------------------------------

MOCK_GNOMAD_TP53_VARIANTS = {
    "data": {
        "gene": {
            "variants": [
                {
                    "variant_id": "17-7675088-C-T",
                    "hgvsp": "p.Arg175His",
                    "consequence": "missense_variant",
                },
                {
                    "variant_id": "17-7675088-C-G",
                    "hgvsp": "p.Arg175Pro",
                    "consequence": "missense_variant",
                },
            ]
        }
    }
}


@respx.mock
async def test_gnomad_find_variant_by_protein_change(http_client):
    respx.post("https://gnomad.broadinstitute.org/api").mock(
        return_value=httpx.Response(200, json=MOCK_GNOMAD_TP53_VARIANTS)
    )
    client = GnomADClient(http_client)
    vid = await client.find_variant_id_by_protein_change("TP53", "p.Arg175His")
    assert vid == "17-7675088-C-T"


@respx.mock
async def test_gnomad_find_variant_missing_returns_none(http_client):
    respx.post("https://gnomad.broadinstitute.org/api").mock(
        return_value=httpx.Response(200, json=MOCK_GNOMAD_TP53_VARIANTS)
    )
    client = GnomADClient(http_client)
    # R175X isn't in the mock → None (cache prevents re-fetch)
    vid = await client.find_variant_id_by_protein_change("TP53", "p.Arg175Cys")
    assert vid is None


@respx.mock
async def test_gnomad_find_variant_network_error_returns_none(http_client):
    respx.post("https://gnomad.broadinstitute.org/api").mock(side_effect=httpx.ConnectError("boom"))
    client = GnomADClient(http_client)
    vid = await client.find_variant_id_by_protein_change("TP53", "p.Arg175His")
    assert vid is None


# ---------------------------------------------------------------------------
# MaveDB per-variant score tests
# ---------------------------------------------------------------------------

MOCK_MAVEDB_BRCA1_SCORES_CSV = (
    "accession,hgvs_nt,hgvs_splice,hgvs_pro,score\n"
    "urn:mavedb:000001#1,NA,NA,p.Arg175His,-2.5\n"
    "urn:mavedb:000001#2,NA,NA,p.Arg175His,-2.3\n"
    "urn:mavedb:000001#3,NA,NA,p.Arg175Cys,-2.8\n"
    "urn:mavedb:000001#4,NA,NA,p.Lys10Ala,0.5\n"
    "urn:mavedb:000001#5,NA,NA,p.Lys10Ala,NA\n"
)


@respx.mock
async def test_mavedb_get_variant_score_match(http_client):
    respx.get("https://api.mavedb.org/api/v1/score-sets/urn:mavedb:000001/scores").mock(
        return_value=httpx.Response(200, text=MOCK_MAVEDB_BRCA1_SCORES_CSV)
    )
    client = MaveDBClient(http_client)
    hits = await client.get_variant_score("urn:mavedb:000001", "p.Arg175His", "Test score set")
    # Two replicates for R175H
    assert len(hits) == 2
    assert all(h.hgvs_pro == "p.Arg175His" for h in hits)
    assert hits[0].score == pytest.approx(-2.5)
    assert hits[0].title == "Test score set"


@respx.mock
async def test_mavedb_get_variant_score_no_match(http_client):
    respx.get("https://api.mavedb.org/api/v1/score-sets/urn:mavedb:000001/scores").mock(
        return_value=httpx.Response(200, text=MOCK_MAVEDB_BRCA1_SCORES_CSV)
    )
    client = MaveDBClient(http_client)
    hits = await client.get_variant_score("urn:mavedb:000001", "p.Val999Ala")
    assert hits == []


@respx.mock
async def test_mavedb_get_variant_score_network_error(http_client):
    respx.get("https://api.mavedb.org/api/v1/score-sets/urn:mavedb:000001/scores").mock(
        side_effect=httpx.ConnectError("boom")
    )
    client = MaveDBClient(http_client)
    hits = await client.get_variant_score("urn:mavedb:000001", "p.Arg175His")
    assert hits == []


@respx.mock
async def test_mavedb_get_variant_score_skips_NA(http_client):
    # Row #5 has score="NA" — must be skipped, not parsed as 0
    respx.get("https://api.mavedb.org/api/v1/score-sets/urn:mavedb:000001/scores").mock(
        return_value=httpx.Response(200, text=MOCK_MAVEDB_BRCA1_SCORES_CSV)
    )
    client = MaveDBClient(http_client)
    hits = await client.get_variant_score("urn:mavedb:000001", "p.Lys10Ala")
    assert len(hits) == 1
    assert hits[0].score == 0.5


# ---------------------------------------------------------------------------
# VariantEffectsClient aggregator tests
# ---------------------------------------------------------------------------


@respx.mock
async def test_variant_effects_aggregator_happy_path(http_client):
    from genesis_bio_mcp.clients.variant_effects import VariantEffectsClient

    # gnomAD gene.variants query
    respx.post("https://gnomad.broadinstitute.org/api").mock(
        return_value=httpx.Response(200, json=MOCK_GNOMAD_TP53_VARIANTS)
    )
    # MyVariant
    respx.get("https://myvariant.info/v1/variant/chr17:g.7675088C>T").mock(
        return_value=httpx.Response(200, json=MOCK_MYVARIANT_TP53_R175H)
    )
    # MaveDB search + scores
    respx.post("https://api.mavedb.org/api/v1/score-sets/search").mock(
        return_value=httpx.Response(
            200,
            json={
                "scoreSets": [
                    {
                        "urn": "urn:mavedb:000001",
                        "title": "TP53 DMS in NSCLC",
                        "numVariants": 5000,
                        "targetGenes": [{"name": "TP53"}],
                        "publishedDate": "2023-01-01",
                    }
                ]
            },
        )
    )
    respx.get("https://api.mavedb.org/api/v1/score-sets/urn:mavedb:000001/scores").mock(
        return_value=httpx.Response(200, text=MOCK_MAVEDB_BRCA1_SCORES_CSV)
    )

    client = VariantEffectsClient(
        gnomad=GnomADClient(http_client),
        myvariant=MyVariantClient(http_client),
        mavedb=MaveDBClient(http_client),
    )
    result = await client.get_effects("TP53", "R175H")

    assert result.gene_symbol == "TP53"
    assert result.canonical_one_letter == "R175H"
    assert result.canonical_hgvs_protein == "p.Arg175His"
    assert result.gnomad_variant_id == "17-7675088-C-T"
    assert result.annotation is not None
    assert result.annotation.clinvar.significance_summary == "Pathogenic"
    assert result.annotation.in_silico.alphamissense_class == "likely_pathogenic"
    assert len(result.dms_scores) == 2  # two R175H rows in the mock CSV
    assert result.notes == []
    md = result.to_markdown()
    assert "R175H" in md
    assert "Pathogenic" in md
    assert "AlphaMissense" in md


@respx.mock
async def test_variant_effects_aggregator_not_in_gnomad(http_client):
    from genesis_bio_mcp.clients.variant_effects import VariantEffectsClient

    respx.post("https://gnomad.broadinstitute.org/api").mock(
        return_value=httpx.Response(200, json={"data": {"gene": {"variants": []}}})
    )
    respx.post("https://api.mavedb.org/api/v1/score-sets/search").mock(
        return_value=httpx.Response(200, json={"scoreSets": []})
    )
    client = VariantEffectsClient(
        gnomad=GnomADClient(http_client),
        myvariant=MyVariantClient(http_client),
        mavedb=MaveDBClient(http_client),
    )
    result = await client.get_effects("TP53", "R175H")

    assert result.gnomad_variant_id is None
    assert result.annotation is None
    assert any("not found in gnomAD" in n for n in result.notes)


async def test_variant_effects_aggregator_invalid_mutation(http_client):
    from genesis_bio_mcp.clients.variant_effects import VariantEffectsClient

    client = VariantEffectsClient(
        gnomad=GnomADClient(http_client),
        myvariant=MyVariantClient(http_client),
        mavedb=MaveDBClient(http_client),
    )
    # parse_protein_change raises ValueError on garbage input
    with pytest.raises(ValueError):
        await client.get_effects("TP53", "not a mutation")
