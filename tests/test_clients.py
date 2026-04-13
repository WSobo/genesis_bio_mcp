"""Unit tests for API client modules using respx to mock httpx."""

import httpx
import pytest
import respx

from genesis_bio_mcp.clients.alphafold import AlphaFoldClient
from genesis_bio_mcp.clients.clinical_trials import ClinicalTrialsClient
from genesis_bio_mcp.clients.depmap import DepMapClient
from genesis_bio_mcp.clients.dgidb import DGIdbClient
from genesis_bio_mcp.clients.gwas import GwasClient
from genesis_bio_mcp.clients.open_targets import OpenTargetsClient
from genesis_bio_mcp.clients.pubchem import PubChemClient
from genesis_bio_mcp.clients.reactome import ReactomeClient
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
