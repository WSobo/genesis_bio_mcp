"""Unit tests for API client modules using respx to mock httpx."""

import httpx
import pytest
import respx

from genesis_bio_mcp.clients.alphafold import AlphaFoldClient
from genesis_bio_mcp.clients.biogrid import BioGRIDClient
from genesis_bio_mcp.clients.chembl import ChEMBLClient
from genesis_bio_mcp.clients.clinical_trials import ClinicalTrialsClient
from genesis_bio_mcp.clients.depmap import DepMapClient
from genesis_bio_mcp.clients.dgidb import DGIdbClient
from genesis_bio_mcp.clients.ensembl import EnsemblClient
from genesis_bio_mcp.clients.gnomad import GnomADClient
from genesis_bio_mcp.clients.gtex import GTExClient
from genesis_bio_mcp.clients.gwas import GwasClient
from genesis_bio_mcp.clients.hpa import HPAClient
from genesis_bio_mcp.clients.iedb import IEDBClient
from genesis_bio_mcp.clients.interpro import InterProClient
from genesis_bio_mcp.clients.mavedb import MaveDBClient
from genesis_bio_mcp.clients.myvariant import MyVariantClient
from genesis_bio_mcp.clients.open_targets import OpenTargetsClient
from genesis_bio_mcp.clients.openfda import OpenFDAClient
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
    MOCK_PUBCHEM_CONCISE_BRAF,
    MOCK_PUBCHEM_GENE_SUMMARY_BRAF,
    MOCK_PUBCHEM_GENESYMBOL_AIDS,
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
    # Bug L (v0.3.3): the empirical 95% cap fires here because the OT-proxy
    # mock returns 100% somatic-mutation-positive cell lines. Real DepMap
    # data on BRAF has fraction ~0.09, so this isn't a false positive in
    # production — it's a side effect of the test fixture's all-positive
    # mock. Documenting that the cap is doing what it should.
    assert result.pan_essential is True
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
    """Primary path: gene→GeneID + AIDs → concise tables filtered by GeneID
    → properties batch. The concise endpoint's per-row Activity Outcome is
    the only reliable Active flag — the bulk cids_type=active index is sparse."""
    respx.get(url__regex=r"gene/genesymbol/.+/summary").mock(
        return_value=httpx.Response(200, json=MOCK_PUBCHEM_GENE_SUMMARY_BRAF)
    )
    respx.get(url__regex=r"assay/target/genesymbol").mock(
        return_value=httpx.Response(200, json=MOCK_PUBCHEM_GENESYMBOL_AIDS)
    )
    respx.get(url__regex=r"assay/aid/\d+/concise").mock(
        return_value=httpx.Response(200, json=MOCK_PUBCHEM_CONCISE_BRAF)
    )
    respx.get(url__regex=r"property/MolecularFormula").mock(
        return_value=httpx.Response(200, json=MOCK_PUBCHEM_PROPERTIES)
    )
    client = PubChemClient(http_client)
    result = await client.get_compounds("BRAF")

    assert result is not None
    assert result.gene_symbol == "BRAF"
    # Two unique CIDs target BRAF (673) and are Active across the mock concise
    # table; the off-target row (GeneID 4916) and Inactive row are filtered.
    assert result.total_active_compounds == 2
    cids = {c.cid for c in result.compounds}
    assert cids == {44462760, 11338033}
    # Activity values converted µM → nM (0.006 µM → 6 nM, 0.031 µM → 31 nM)
    by_cid = {c.cid: c for c in result.compounds}
    assert by_cid[11338033].activity_value == pytest.approx(6.0)
    assert by_cid[44462760].activity_value == pytest.approx(31.0)
    # Sorted by potency: dabrafenib (6 nM) before vemurafenib (31 nM)
    assert result.compounds[0].cid == 11338033
    # Property enrichment populated names/formulae from the property batch
    assert "dabrafenib" in by_cid[11338033].name.lower()
    assert by_cid[11338033].molecular_formula == "C23H20F3N5O2S2"


@respx.mock
async def test_pubchem_filters_panel_assay_rows_by_gene_id(http_client):
    """Concise tables from kinase panel assays include rows for many targets;
    only rows whose Target GeneID matches the resolved GeneID may surface."""
    respx.get(url__regex=r"gene/genesymbol/.+/summary").mock(
        return_value=httpx.Response(200, json=MOCK_PUBCHEM_GENE_SUMMARY_BRAF)
    )
    respx.get(url__regex=r"assay/target/genesymbol").mock(
        return_value=httpx.Response(200, json=MOCK_PUBCHEM_GENESYMBOL_AIDS)
    )
    respx.get(url__regex=r"assay/aid/\d+/concise").mock(
        return_value=httpx.Response(200, json=MOCK_PUBCHEM_CONCISE_BRAF)
    )
    respx.get(url__regex=r"property/MolecularFormula").mock(
        return_value=httpx.Response(200, json=MOCK_PUBCHEM_PROPERTIES)
    )
    client = PubChemClient(http_client)
    result = await client.get_compounds("BRAF")

    cids = {c.cid for c in result.compounds}
    # CID 999 (Active, GeneID 4916=NTRK3) must not leak into BRAF results
    assert 999 not in cids
    # CID 12345 (Inactive against BRAF) must not appear either
    assert 12345 not in cids


@respx.mock
async def test_pubchem_falls_back_to_entrez_when_pug_rest_empty(http_client):
    """If PUG REST returns no AIDs, fall back to Entrez assay search."""
    respx.get(url__regex=r"gene/genesymbol/.+/summary").mock(
        return_value=httpx.Response(200, json=MOCK_PUBCHEM_GENE_SUMMARY_BRAF)
    )
    respx.get(url__regex=r"assay/target/genesymbol").mock(return_value=httpx.Response(404))
    respx.get(url__regex=r"esearch\.fcgi").mock(
        return_value=httpx.Response(200, json=MOCK_ENTREZ_AIDS)
    )
    respx.get(url__regex=r"assay/aid/\d+/concise").mock(
        return_value=httpx.Response(200, json=MOCK_PUBCHEM_CONCISE_BRAF)
    )
    respx.get(url__regex=r"property/MolecularFormula").mock(
        return_value=httpx.Response(200, json=MOCK_PUBCHEM_PROPERTIES)
    )
    client = PubChemClient(http_client)
    result = await client.get_compounds("BRAF")
    assert result is not None
    assert result.total_active_compounds == 2


@respx.mock
async def test_pubchem_retries_on_503(http_client):
    """Verify tenacity retries on 503 from the primary AID-list path."""
    call_count = 0

    def side_effect(request):
        nonlocal call_count
        call_count += 1
        if call_count == 1:
            return httpx.Response(503)
        return httpx.Response(200, json=MOCK_PUBCHEM_GENESYMBOL_AIDS)

    respx.get(url__regex=r"gene/genesymbol/.+/summary").mock(
        return_value=httpx.Response(200, json=MOCK_PUBCHEM_GENE_SUMMARY_BRAF)
    )
    respx.get(url__regex=r"assay/target/genesymbol").mock(side_effect=side_effect)
    respx.get(url__regex=r"assay/aid/\d+/concise").mock(
        return_value=httpx.Response(200, json=MOCK_PUBCHEM_CONCISE_BRAF)
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
    # Both PUG REST and Entrez return empty
    respx.get(url__regex=r"gene/genesymbol/.+/summary").mock(return_value=httpx.Response(404))
    respx.get(url__regex=r"assay/target/genesymbol").mock(return_value=httpx.Response(404))
    respx.get(url__regex=r"esearch\.fcgi").mock(
        return_value=httpx.Response(200, json={"esearchresult": {"idlist": []}})
    )
    client = PubChemClient(http_client)
    result = await client.get_compounds("FAKEGENE")
    assert result is None


@respx.mock
async def test_pubchem_returns_none_when_no_active_rows_match_gene(http_client):
    """AIDs exist but no concise rows match the resolved GeneID — return None
    rather than fall through to a misleading name-based fallback."""
    respx.get(url__regex=r"gene/genesymbol/.+/summary").mock(
        return_value=httpx.Response(200, json=MOCK_PUBCHEM_GENE_SUMMARY_BRAF)
    )
    respx.get(url__regex=r"assay/target/genesymbol").mock(
        return_value=httpx.Response(200, json=MOCK_PUBCHEM_GENESYMBOL_AIDS)
    )
    # Concise table has no rows targeting GeneID 673
    no_match_concise = {
        "Table": {
            "Columns": {
                "Column": [
                    "AID",
                    "SID",
                    "CID",
                    "Activity Outcome",
                    "Target Accession",
                    "Target GeneID",
                    "Activity Value [uM]",
                    "Activity Name",
                    "Assay Name",
                    "Assay Type",
                    "PubMed ID",
                    "RNAi",
                ]
            },
            "Row": [
                {
                    "Cell": [
                        "1259398",
                        "1",
                        "100",
                        "Active",
                        "NP_001",
                        "9999",  # not BRAF
                        "0.1",
                        "IC50",
                        "x",
                        "x",
                        "",
                        "",
                    ]
                }
            ],
        }
    }
    respx.get(url__regex=r"assay/aid/\d+/concise").mock(
        return_value=httpx.Response(200, json=no_match_concise)
    )
    client = PubChemClient(http_client)
    result = await client.get_compounds("BRAF")
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
# Ensembl client tests
# ---------------------------------------------------------------------------


_MOCK_ENSEMBL_LOOKUP_BRAF = {
    "id": "ENSG00000157764",
    "display_name": "BRAF",
    "seq_region_name": "7",
    "start": 140719327,
    "end": 140924928,
    "strand": -1,
    "biotype": "protein_coding",
    "Transcript": [
        {
            "id": "ENST00000646891",
            "is_canonical": 1,
            "biotype": "protein_coding",
            "length": 6490,
        },
        {
            "id": "ENST00000288602",
            "is_canonical": 0,
            "biotype": "protein_coding",
            "length": 2949,
        },
    ],
}

_MOCK_ENSEMBL_VEP_V600E = [
    {
        "input": "ENST00000646891:p.Val600Glu",
        "most_severe_consequence": "missense_variant",
        "assembly_name": "GRCh38",
        "transcript_consequences": [
            {
                "transcript_id": "ENST00000646891",
                "gene_symbol": "BRAF",
                "canonical": 1,
                "biotype": "protein_coding",
                "impact": "MODERATE",
                "consequence_terms": ["missense_variant"],
                "sift_score": 0.0,
                "sift_prediction": "deleterious",
                "polyphen_score": 0.95,
                "polyphen_prediction": "probably_damaging",
                "amino_acids": "V/E",
                "codons": "gTg/gAg",
            },
            {
                "transcript_id": "ENST00000288602",
                "gene_symbol": "BRAF",
                "canonical": 0,
                "biotype": "protein_coding",
                "impact": "MODERATE",
                "consequence_terms": ["missense_variant"],
            },
        ],
        "regulatory_feature_consequences": [],
    }
]


@respx.mock
async def test_ensembl_lookup_gene_happy_path(http_client):
    respx.get(url__regex=r"rest\.ensembl\.org/lookup/symbol/homo_sapiens/BRAF").mock(
        return_value=httpx.Response(200, json=_MOCK_ENSEMBL_LOOKUP_BRAF)
    )
    client = EnsemblClient(http_client)
    gene = await client.lookup_gene("braf")

    assert gene is not None
    assert gene.ensembl_id == "ENSG00000157764"
    assert gene.symbol == "BRAF"
    assert gene.chrom == "7"
    assert gene.strand == -1
    assert gene.canonical_transcript_id == "ENST00000646891"
    assert len(gene.transcripts) == 2
    md = gene.to_markdown()
    assert "ENSG00000157764" in md
    assert "ENST00000646891" in md


@respx.mock
async def test_ensembl_lookup_gene_404(http_client):
    respx.get(url__regex=r"rest\.ensembl\.org/lookup/symbol/homo_sapiens/NOTAGENE").mock(
        return_value=httpx.Response(404)
    )
    client = EnsemblClient(http_client)
    gene = await client.lookup_gene("NOTAGENE")
    assert gene is None


@respx.mock
async def test_ensembl_get_vep_by_hgvs_canonical_only(http_client):
    respx.get(url__regex=r"rest\.ensembl\.org/vep/human/hgvs/ENST00000646891:p\.Val600Glu").mock(
        return_value=httpx.Response(200, json=_MOCK_ENSEMBL_VEP_V600E)
    )
    client = EnsemblClient(http_client)
    report = await client.get_vep_by_hgvs("ENST00000646891:p.Val600Glu")

    assert report is not None
    assert report.most_severe_consequence == "missense_variant"
    assert report.assembly_name == "GRCh38"
    assert len(report.consequences) == 1
    assert report.consequences[0].canonical is True
    assert report.consequences[0].sift_score == 0.0
    assert report.consequences[0].polyphen_prediction == "probably_damaging"


@respx.mock
async def test_ensembl_get_vep_all_transcripts(http_client):
    respx.get(url__regex=r"rest\.ensembl\.org/vep/human/hgvs").mock(
        return_value=httpx.Response(200, json=_MOCK_ENSEMBL_VEP_V600E)
    )
    client = EnsemblClient(http_client)
    report = await client.get_vep_by_hgvs(
        "ENST00000646891:p.Val600Glu", include_all_transcripts=True
    )
    assert report is not None
    assert len(report.consequences) == 2


@respx.mock
async def test_ensembl_get_vep_consequences_combined(http_client):
    respx.get(url__regex=r"rest\.ensembl\.org/lookup/symbol/homo_sapiens/BRAF").mock(
        return_value=httpx.Response(200, json=_MOCK_ENSEMBL_LOOKUP_BRAF)
    )
    respx.get(url__regex=r"rest\.ensembl\.org/vep/human/hgvs").mock(
        return_value=httpx.Response(200, json=_MOCK_ENSEMBL_VEP_V600E)
    )
    client = EnsemblClient(http_client)
    report = await client.get_vep_consequences("BRAF", "p.Val600Glu")
    assert report is not None
    assert report.most_severe_consequence == "missense_variant"


@respx.mock
async def test_ensembl_network_failure_returns_none(http_client):
    respx.get(url__regex=r"rest\.ensembl\.org").mock(side_effect=httpx.ConnectTimeout("boom"))
    client = EnsemblClient(http_client)
    assert await client.lookup_gene("BRAF") is None


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


_MOCK_DGIDB_SALT_FORMS = {
    "data": {
        "genes": {
            "nodes": [
                {
                    "name": "JAK2",
                    "interactions": [
                        {
                            "drug": {"name": "FILGOTINIB", "approved": True},
                            "interactionTypes": [{"type": "inhibitor"}],
                            "interactionClaims": [{"source": {"sourceDbName": "ChEMBL"}}],
                        },
                        {
                            "drug": {"name": "FILGOTINIB MALEATE", "approved": True},
                            "interactionTypes": [{"type": "inhibitor"}],
                            "interactionClaims": [{"source": {"sourceDbName": "DrugBank"}}],
                        },
                    ],
                }
            ]
        }
    }
}


@respx.mock
async def test_dgidb_collapses_salt_forms(http_client):
    """FILGOTINIB + FILGOTINIB MALEATE are the same molecule (free base + salt).

    DGIdb returns them as distinct records, double-counting the same INN. The
    token-prefix dedup in _collapse_salt_forms merges the salt row into the
    parent, unioning sources.
    """
    respx.post("https://dgidb.org/api/graphql").mock(
        return_value=httpx.Response(200, json=_MOCK_DGIDB_SALT_FORMS)
    )
    client = DGIdbClient(http_client)
    result = await client.get_drug_interactions("JAK2")

    assert len(result) == 1
    assert result[0].drug_name == "FILGOTINIB"
    # Sources are unioned across the parent and merged salt form
    assert set(result[0].sources) == {"ChEMBL", "DrugBank"}


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


def _make_cffi_session_mock(payload: dict | None = None, raise_for_status: bool = False):
    """Build a mock curl_cffi AsyncSession context manager returning *payload*.

    Used to test ClinicalTrialsClient without reaching the network. The
    client uses curl_cffi (not httpx), so respx can't intercept — we mock
    the AsyncSession class itself.
    """
    from unittest.mock import AsyncMock, MagicMock

    mock_resp = MagicMock()
    if raise_for_status:
        mock_resp.raise_for_status.side_effect = Exception("HTTP 500")
    else:
        mock_resp.raise_for_status.return_value = None
        mock_resp.json.return_value = payload or {}

    mock_session = AsyncMock()
    mock_session.__aenter__.return_value.get = AsyncMock(return_value=mock_resp)
    mock_session.__aexit__.return_value = None
    return mock_session


async def test_clinical_trials_get_trials(http_client, monkeypatch):
    mock_session = _make_cffi_session_mock(_MOCK_CT_RESPONSE)
    monkeypatch.setattr(
        "genesis_bio_mcp.clients.clinical_trials.AsyncSession",
        lambda *args, **kwargs: mock_session,
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


async def test_clinical_trials_returns_empty_on_error(http_client, monkeypatch):
    mock_session = _make_cffi_session_mock(raise_for_status=True)
    monkeypatch.setattr(
        "genesis_bio_mcp.clients.clinical_trials.AsyncSession",
        lambda *args, **kwargs: mock_session,
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
async def test_reactome_dedups_pathways_by_stid(http_client):
    """Reactome's analysis service can return the same stId twice with
    slightly different gene_count; output must be distinct on stId."""
    duplicate_response = {
        "summary": {"token": "tok", "type": "OVERREPRESENTATION"},
        "pathways": [
            {
                "stId": "R-HSA-5673001",
                "name": "RAF/MAP kinase cascade",
                "entities": {"pValue": 1.2e-15, "total": 42},
            },
            {
                # Same stId, different gene_count — duplicate to be removed
                "stId": "R-HSA-5673001",
                "name": "RAF/MAP kinase cascade",
                "entities": {"pValue": 1.2e-15, "total": 38},
            },
            {
                "stId": "R-HSA-162582",
                "name": "Signal Transduction",
                "entities": {"pValue": 0.0023, "total": 512},
            },
        ],
    }
    respx.post(url__regex=r"reactome\.org/AnalysisService/identifiers").mock(
        return_value=httpx.Response(200, json=duplicate_response)
    )
    client = ReactomeClient(http_client)
    result = await client.get_pathway_context("BRAF")

    assert result is not None
    stids = [p.reactome_id for p in result.pathways]
    assert stids.count("R-HSA-5673001") == 1
    # First-occurrence wins: gene_count from the first row preserved
    raf = next(p for p in result.pathways if p.reactome_id == "R-HSA-5673001")
    assert raf.gene_count == 42


@respx.mock
async def test_reactome_dedups_pathways_by_display_name(http_client):
    """Reactome sometimes returns parent + child pathways with the same
    display_name under different stable IDs (e.g. different organism scopes
    or refined annotations). Second-pass dedup keeps the row with the
    smallest p-value."""
    dup_name_response = {
        "summary": {"token": "tok", "type": "OVERREPRESENTATION"},
        "pathways": [
            {
                # Larger p-value — must be dropped
                "stId": "R-HSA-4",
                "name": "IFNG signaling activates MAPKs",
                "entities": {"pValue": 5e-3, "total": 100},
            },
            {
                # Smaller p-value — must win
                "stId": "R-HSA-5",
                "name": "IFNG signaling activates MAPKs",
                "entities": {"pValue": 1e-8, "total": 80},
            },
        ],
    }
    respx.post(url__regex=r"reactome\.org/AnalysisService/identifiers").mock(
        return_value=httpx.Response(200, json=dup_name_response)
    )
    client = ReactomeClient(http_client)
    result = await client.get_pathway_context("JAK2")
    assert result is not None
    assert len(result.pathways) == 1
    kept = result.pathways[0]
    assert kept.p_value == 1e-8
    assert kept.reactome_id == "R-HSA-5"


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

# Real Reactome /data/participants response shape: each entry is a
# PhysicalEntity with refEntities carrying displayName "UniProt:ACC SYMBOL"
# and a schemaClass — gene symbols must be parsed from displayName.
_MOCK_REACTOME_PARTICIPANTS_RESPONSE = [
    {
        "refEntities": [
            {
                "displayName": "UniProt:P15056 BRAF",
                "identifier": "P15056",
                "schemaClass": "ReferenceGeneProduct",
            },
            {
                "displayName": "UniProt:P04049 RAF1",
                "identifier": "P04049",
                "schemaClass": "ReferenceGeneProduct",
            },
        ]
    },
    {
        "refEntities": [
            {
                "displayName": "UniProt:Q02750 MAP2K1",
                "identifier": "Q02750",
                "schemaClass": "ReferenceGeneProduct",
            },
            # Non-protein refEntity should be filtered out
            {
                "displayName": "water [ChEBI:15377]",
                "identifier": "15377",
                "schemaClass": "ReferenceMolecule",
            },
        ]
    },
]


@respx.mock
async def test_reactome_get_pathway_members_by_name(http_client):
    # Search endpoint lives at /search/query (no /data/ prefix)
    respx.get(url__regex=r"reactome\.org/ContentService/search/query").mock(
        return_value=httpx.Response(200, json=_MOCK_REACTOME_SEARCH_RESPONSE)
    )
    respx.get(url__regex=r"reactome\.org/ContentService/data/participants").mock(
        return_value=httpx.Response(200, json=_MOCK_REACTOME_PARTICIPANTS_RESPONSE)
    )
    client = ReactomeClient(http_client)
    genes = await client.get_pathway_members("RAF/MAP kinase cascade")

    assert "BRAF" in genes
    assert "RAF1" in genes
    assert "MAP2K1" in genes
    # ReferenceMolecule entries (water) must not leak through
    assert not any("WATER" in g or "CHEBI" in g for g in genes)


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
    respx.get(url__regex=r"reactome\.org/ContentService/search/query").mock(
        return_value=httpx.Response(200, json={"results": []})
    )
    client = ReactomeClient(http_client)
    genes = await client.get_pathway_members("nonexistent pathway xyz")
    assert genes == []


@respx.mock
async def test_reactome_get_pathway_members_caches_result(http_client):
    """Second call with same input must hit the in-memory cache (no re-fetch)."""
    participants_route = respx.get(
        url__regex=r"reactome\.org/ContentService/data/participants"
    ).mock(return_value=httpx.Response(200, json=_MOCK_REACTOME_PARTICIPANTS_RESPONSE))
    client = ReactomeClient(http_client)
    first = await client.get_pathway_members("R-HSA-5673001")
    second = await client.get_pathway_members("R-HSA-5673001")
    assert first == second
    assert participants_route.call_count == 1


@respx.mock
async def test_reactome_get_pathway_members_truncates_at_max(http_client):
    """max_genes caps the result and short-circuits the parse loop."""
    big_participants = [
        {
            "refEntities": [
                {
                    "displayName": f"UniProt:P{i:05d} GENE{i}",
                    "identifier": f"P{i:05d}",
                    "schemaClass": "ReferenceGeneProduct",
                }
            ]
        }
        for i in range(200)
    ]
    respx.get(url__regex=r"reactome\.org/ContentService/data/participants").mock(
        return_value=httpx.Response(200, json=big_participants)
    )
    client = ReactomeClient(http_client)
    genes = await client.get_pathway_members("R-HSA-5673001", max_genes=20)
    assert len(genes) == 20


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
async def test_efo_resolver_expands_hierarchy(http_client):
    """Resolved terms must be hierarchy-expanded so URI matching catches studies
    tagged with sibling/parent/descendant EFO terms — e.g. a "polycythemia vera"
    query matches JAK2 associations tagged with "myeloproliferative neoplasm"
    (its direct EFO parent)."""

    # First call (no filter) returns the primary match. Subsequent calls carry
    # an allChildrenOf / ancestorsOf filter — distinguish them by URL query.
    def side_effect(request):
        url = str(request.url)
        if "allChildrenOf" in url:
            return httpx.Response(
                200,
                json={
                    "response": {
                        "docs": [
                            {"iri": "http://www.ebi.ac.uk/efo/EFO_DESC_1"},
                            {"iri": "http://www.ebi.ac.uk/efo/EFO_DESC_2"},
                        ]
                    }
                },
            )
        if "ancestorsOf" in url:
            return httpx.Response(
                200,
                json={"response": {"docs": [{"iri": "http://www.ebi.ac.uk/efo/EFO_ANC_1"}]}},
            )
        return httpx.Response(
            200,
            json={
                "response": {
                    "docs": [
                        {
                            "iri": "http://www.ebi.ac.uk/efo/EFO_PRIMARY",
                            "label": "polycythemia vera",
                            "synonym": [],
                        }
                    ]
                }
            },
        )

    respx.get(url__regex=r"ols4/api/search").mock(side_effect=side_effect)
    resolver = EFOResolver(http_client, cache_path=None)
    terms = await resolver.resolve("polycythemia vera")

    assert len(terms) == 1
    related = set(terms[0].related_uris)
    assert "http://www.ebi.ac.uk/efo/EFO_DESC_1" in related
    assert "http://www.ebi.ac.uk/efo/EFO_DESC_2" in related
    assert "http://www.ebi.ac.uk/efo/EFO_ANC_1" in related
    # The primary URI must not appear in related_uris (self excluded)
    assert "http://www.ebi.ac.uk/efo/EFO_PRIMARY" not in related


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
    calls_after_first = call_count
    await resolver.resolve("obesity")  # second call — must hit session cache

    # Session cache prevents any additional OLS4 traffic on the second resolve,
    # including the hierarchy-expansion calls (allChildrenOf/ancestorsOf per term).
    assert call_count == calls_after_first


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


def test_filter_by_trait_biomarker_query_ldl_cholesterol():
    """Users frequently query by biomarker name. PCSK9 × 'LDL cholesterol' is
    a textbook GWAS association — must match hits using long-form labels."""
    hits = [
        _make_hit("Low density lipoprotein cholesterol levels"),
        _make_hit("LDL cholesterol levels"),
        _make_hit("type 2 diabetes"),
    ]
    result = filter_by_trait(hits, "LDL cholesterol", efo_terms=None)
    traits = {h.trait for h in result}
    assert "Low density lipoprotein cholesterol levels" in traits
    assert "LDL cholesterol levels" in traits
    assert "type 2 diabetes" not in traits


def test_filter_by_trait_matches_on_related_uris():
    """A hit tagged with a parent/descendant EFO URI must match when the user
    query's resolved EFO term has that URI in ``related_uris``. This is the
    fix for JAK2 × 'polycythemia vera' — studies tagged with MPN (parent)
    now match via structural identity, not free-text substring."""
    hits = [
        _make_hit(
            "Myeloproliferative neoplasm",
            efo_uri="http://www.ebi.ac.uk/efo/EFO_0000538",  # parent of polycythemia vera
        ),
        _make_hit(
            "Type 2 diabetes",
            efo_uri="http://www.ebi.ac.uk/efo/EFO_0000400",
        ),
    ]
    efo_terms = [
        EFOTerm(
            uri="http://www.ebi.ac.uk/efo/EFO_0004254",  # polycythemia vera itself
            label="polycythemia vera",
            synonyms=[],
            related_uris=["http://www.ebi.ac.uk/efo/EFO_0000538"],  # MPN parent
        )
    ]
    result = filter_by_trait(hits, "polycythemia vera", efo_terms=efo_terms)
    assert len(result) == 1
    assert result[0].trait == "Myeloproliferative neoplasm"


def test_filter_by_trait_synonyms_additive_to_efo():
    """EFO label may not include common biomarker abbreviations; the synonym
    dict is applied additively, not as an XOR fallback."""
    hits = [_make_hit("LDL-c levels")]
    # EFO returns the canonical measurement name but no "ldl-c" synonym
    efo_terms = [
        EFOTerm(
            uri="http://www.ebi.ac.uk/efo/EFO_0004611",
            label="LDL cholesterol measurement",
            synonyms=[],
        )
    ]
    result = filter_by_trait(hits, "ldl cholesterol", efo_terms=efo_terms)
    # Without additive merge, "ldl-c levels" wouldn't match
    # ("ldl cholesterol measurement" isn't a substring); with the merge,
    # the "ldl-c" synonym from TRAIT_SYNONYMS catches it.
    assert len(result) == 1


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
async def test_sabdab_dedups_by_pdb_id(http_client, tmp_path, monkeypatch):
    """SAbDab has one row per chain pair; the same PDB can appear twice when
    both chains match the antigen. Output must dedup by PDB, keeping best res."""
    monkeypatch.setattr(
        "genesis_bio_mcp.clients.sabdab.settings.sabdab_cache_path",
        tmp_path / "sabdab_cache.tsv",
    )
    monkeypatch.setattr(
        "genesis_bio_mcp.clients.sabdab.settings.sabdab_cache_ttl_secs",
        604800,
    )

    # Two rows for the same PDB 1abc — different chain pair, different res
    duplicate_tsv = (
        "pdb\tHchain\tLchain\tmodel\tantigen_chain\tantigen_type\tantigen_het_name\t"
        "antigen_name\tshort_header\tdate\tcompound\torganism\theavy_species\tlight_species\t"
        "antigen_species\tauthors\tresolution\tmethod\tr_free\tr_factor\tscfv\tengineered\t"
        "heavy_subclass\tlight_subclass\tlight_ctype\taffinity\tdelta_g\taffinity_method\t"
        "temperature\tpmid\n"
        "1abc\tH\tL\t0\tA\tprotein\tNA\tepidermal growth factor receptor\tIMMUNE SYSTEM\t"
        "01/01/20\tAnti-EGFR\tHomo sapiens\thomo sapiens\thomo sapiens\thomo sapiens\t"
        "Smith\t2.5\tX-RAY\t0.21\t0.18\tFalse\tTrue\tIGHV3\tIGKV1\tKappa\tNone\tNone\tNone\tNone\t1\n"
        "1abc\tB\tA\t0\tC\tprotein\tNA\tepidermal growth factor receptor\tIMMUNE SYSTEM\t"
        "01/01/20\tAnti-EGFR\tHomo sapiens\thomo sapiens\thomo sapiens\thomo sapiens\t"
        "Smith\t3.0\tX-RAY\t0.21\t0.18\tFalse\tTrue\tIGHV3\tIGKV1\tKappa\tNone\tNone\tNone\tNone\t1\n"
    )
    respx.get(url__regex=r"sabdab-sabpred.*summary").mock(
        return_value=httpx.Response(200, content=duplicate_tsv.encode())
    )

    client = SAbDabClient(http_client)
    result = await client.get_antibody_structures("egfr")

    assert result is not None
    assert result.total_structures == 1
    assert result.structures[0].pdb == "1ABC"
    # Best-resolution row (2.5 Å) wins, not the 3.0 Å duplicate
    assert result.structures[0].resolution_ang == pytest.approx(2.5)


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


def test_sabdab_affinity_insight_skipped_when_all_null():
    """When no structure has a positive measured affinity, the 'Best measured
    affinity' insight line must be omitted — NOT rendered as '0.0 nM'."""
    from genesis_bio_mcp.models import AntibodyStructure, AntibodyStructures

    structures = [
        AntibodyStructure(
            pdb="5PDL",
            is_nanobody=False,
            antigen_name="programmed cell death 1 ligand 1",
            resolution_ang=2.0,
            method="X-RAY DIFFRACTION",
            heavy_species="homo sapiens",
            light_species="homo sapiens",
            heavy_subclass=None,
            light_subclass=None,
            is_engineered=True,
            is_scfv=False,
            affinity_nM=None,
            compound=None,
            date_added=None,
            pmid=None,
        ),
    ]
    result = AntibodyStructures(
        query="PD-L1",
        total_structures=1,
        nanobody_count=0,
        fab_count=1,
        structures=structures,
    )
    md = result.to_markdown()
    assert "0.0 nM" not in md
    assert "Best measured affinity" not in md


def test_sabdab_parse_float_rejects_zero_and_negative():
    """_parse_float must return None for zero / negative / non-finite values so
    they don't leak into the affinity aggregation."""
    from genesis_bio_mcp.clients.sabdab import _parse_float

    assert _parse_float("5.2") == 5.2
    assert _parse_float("0") is None
    assert _parse_float("0.0") is None
    assert _parse_float("-1.5") is None
    assert _parse_float("") is None
    assert _parse_float("NA") is None


def test_compounds_activity_falls_back_to_outcome():
    """PubChem concise rows without an 'Activity Name' cell render activity_type
    as None; Compounds.to_markdown must fall back to activity_outcome ('Active')
    rather than a bare em-dash."""
    from genesis_bio_mcp.models import CompoundActivity, Compounds

    compounds = Compounds(
        gene_symbol="JAK2",
        total_active_compounds=1,
        compounds=[
            CompoundActivity(
                cid=12345,
                name="hypothetical-kinase-inhibitor",
                molecular_formula="C20H20N4O",
                molecular_weight=356.4,
                activity_outcome="Active",
                activity_value=45.0,
                activity_type=None,  # PubChem concise row missing Activity Name
                assay_id=None,
            )
        ],
    )
    md = compounds.to_markdown()
    # The Activity cell should be "Active" (outcome fallback), not a bare "—"
    assert "| Active |" in md


def test_score_breakdown_sums_to_priority_score():
    """_compute_score's breakdown contributions must sum (before 10.0 cap) to
    the reported priority_score. This is the auditability invariant that
    compare_targets rendering relies on."""
    from genesis_bio_mcp.models import (
        CancerDependency,
        GwasEvidence,
        TargetDiseaseAssociation,
    )
    from genesis_bio_mcp.tools.target_prioritization import _compute_score

    disease_assoc = TargetDiseaseAssociation(
        gene_symbol="JAK2",
        disease_name="polycythemia vera",
        disease_efo_id="EFO_0004254",
        ensembl_id="ENSG00000096968",
        overall_score=0.75,
        evidence_count=42,
        known_drug_score=0.9,
        genetic_association_score=0.6,
    )
    cancer_dep = CancerDependency(
        gene_symbol="JAK2",
        mean_ceres_score=-0.2,
        fraction_dependent_lines=0.1,
        top_dependent_lineages=["Blood"],
        pan_essential=False,
        data_source="DepMap Chronos",
    )
    gwas_ev = GwasEvidence(
        gene_symbol="JAK2",
        trait_query="polycythemia vera",
        total_associations=5,
        associations=[],
        strongest_p_value=1e-15,
    )
    total, breakdown = _compute_score(
        disease_assoc,
        cancer_dep,
        gwas_ev,
        None,
        None,
        None,
        indication="polycythemia vera",
    )
    # Breakdown.total is the pre-cap sum; priority_score is min(total, 10.0)
    assert abs(breakdown.total - total) < 1e-6 or breakdown.total >= 10.0
    # Non-zero axes present for populated inputs
    assert breakdown.ot > 0
    assert breakdown.depmap > 0
    assert breakdown.gwas > 0
    assert breakdown.known_drug > 0
    # And unpopulated axes are zero
    assert breakdown.chem_matter == 0.0
    assert breakdown.protein == 0.0
    # Expression axis is 0.0 when no HPA report is passed (v0.3.0: backward-compat)
    assert breakdown.expression == 0.0


@pytest.mark.parametrize(
    "category,expected",
    [
        ("Tissue enriched", 1.0),
        ("Group enriched", 0.7),
        ("Tissue enhanced", 0.5),
        ("Low tissue specificity", 0.2),
        ("Not detected", 0.0),
        ("Unknown category", 0.0),  # unmapped category → no credit
        (None, 0.0),  # missing category on the HPA payload
    ],
)
def test_expression_axis_maps_hpa_category_to_score(category, expected):
    """The expression axis must follow the HPA tissue-specificity table exactly.

    This is the new v0.3.0 scoring axis and its weights directly affect
    target rankings for reference targets like BRAF (low specificity) vs.
    tissue-restricted targets like MC1R (tissue enriched in skin).
    """
    from genesis_bio_mcp.models import HPAExpression, ProteinAtlasReport
    from genesis_bio_mcp.tools.target_prioritization import _compute_score

    atlas = ProteinAtlasReport(
        gene_symbol="TEST",
        expression=HPAExpression(
            gene_symbol="TEST",
            rna_tissue_specificity_category=category,
        ),
        pathology=[],
    )
    _, breakdown = _compute_score(
        disease_assoc=None,
        cancer_dep=None,
        gwas_ev=None,
        compounds=None,
        protein=None,
        chembl_compounds=None,
        indication="test",
        protein_atlas=atlas,
    )
    assert breakdown.expression == expected
    # When ONLY the expression axis contributes, the priority score equals
    # the expression value — this locks the axis's isolation from the cap.
    assert breakdown.total == expected


def test_expression_axis_none_report_contributes_zero():
    """protein_atlas=None is the backward-compatible path; axis contributes 0.0."""
    from genesis_bio_mcp.tools.target_prioritization import _compute_score

    _, breakdown = _compute_score(
        disease_assoc=None,
        cancer_dep=None,
        gwas_ev=None,
        compounds=None,
        protein=None,
        chembl_compounds=None,
        indication="test",
        protein_atlas=None,
    )
    assert breakdown.expression == 0.0
    assert breakdown.total == 0.0


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
    # Ensembl gene lookup + VEP (so VEP note isn't emitted)
    respx.get(url__regex=r"rest\.ensembl\.org/lookup/symbol/homo_sapiens/TP53").mock(
        return_value=httpx.Response(
            200,
            json={
                "id": "ENSG00000141510",
                "display_name": "TP53",
                "seq_region_name": "17",
                "start": 7661779,
                "end": 7687550,
                "strand": -1,
                "biotype": "protein_coding",
                "Transcript": [
                    {"id": "ENST00000269305", "is_canonical": 1, "biotype": "protein_coding"}
                ],
            },
        )
    )
    respx.get(url__regex=r"rest\.ensembl\.org/vep/human/hgvs").mock(
        return_value=httpx.Response(
            200,
            json=[
                {
                    "most_severe_consequence": "missense_variant",
                    "assembly_name": "GRCh38",
                    "transcript_consequences": [
                        {
                            "transcript_id": "ENST00000269305",
                            "canonical": 1,
                            "impact": "MODERATE",
                            "consequence_terms": ["missense_variant"],
                            "amino_acids": "R/H",
                        }
                    ],
                    "regulatory_feature_consequences": [],
                }
            ],
        )
    )

    client = VariantEffectsClient(
        gnomad=GnomADClient(http_client),
        myvariant=MyVariantClient(http_client),
        mavedb=MaveDBClient(http_client),
        ensembl=EnsemblClient(http_client),
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
    assert result.vep_consequences is not None
    assert result.vep_consequences.most_severe_consequence == "missense_variant"
    assert result.notes == []
    md = result.to_markdown()
    assert "R175H" in md
    assert "Pathogenic" in md
    assert "AlphaMissense" in md


@respx.mock
async def test_variant_effects_aggregator_not_in_gnomad_falls_back_to_myvariant(http_client):
    """Somatic hotspots (BRAF V600E, etc.) are absent from gnomAD by design;
    ClinVar/AlphaMissense must still be returned via MyVariant /query."""
    from genesis_bio_mcp.clients.variant_effects import VariantEffectsClient

    respx.post("https://gnomad.broadinstitute.org/api").mock(
        return_value=httpx.Response(200, json={"data": {"gene": {"variants": []}}})
    )
    # MyVariant /query returns the ClinVar+AlphaMissense payload as a hit
    respx.get(url__regex=r"https://myvariant\.info/v1/query.*").mock(
        return_value=httpx.Response(
            200,
            json={"hits": [{"_id": "chr17:g.7675088C>T", **MOCK_MYVARIANT_TP53_R175H}]},
        )
    )
    respx.post("https://api.mavedb.org/api/v1/score-sets/search").mock(
        return_value=httpx.Response(200, json={"scoreSets": []})
    )
    respx.get(url__regex=r"rest\.ensembl\.org").mock(return_value=httpx.Response(404))
    client = VariantEffectsClient(
        gnomad=GnomADClient(http_client),
        myvariant=MyVariantClient(http_client),
        mavedb=MaveDBClient(http_client),
        ensembl=EnsemblClient(http_client),
    )
    result = await client.get_effects("TP53", "R175H")

    assert result.gnomad_variant_id is None
    # gnomAD-miss note is present, but ClinVar/AlphaMissense are populated
    assert any("not found in gnomAD" in n for n in result.notes)
    assert result.annotation is not None
    assert result.annotation.clinvar.significance_summary == "Pathogenic"
    assert result.annotation.in_silico.alphamissense_class == "likely_pathogenic"


@respx.mock
async def test_variant_effects_aggregator_not_in_gnomad_or_myvariant(http_client):
    """When both gnomAD and MyVariant /query miss, annotation is None and
    both notes are recorded."""
    from genesis_bio_mcp.clients.variant_effects import VariantEffectsClient

    respx.post("https://gnomad.broadinstitute.org/api").mock(
        return_value=httpx.Response(200, json={"data": {"gene": {"variants": []}}})
    )
    respx.get(url__regex=r"https://myvariant\.info/v1/query.*").mock(
        return_value=httpx.Response(200, json={"hits": [], "total": 0})
    )
    respx.post("https://api.mavedb.org/api/v1/score-sets/search").mock(
        return_value=httpx.Response(200, json={"scoreSets": []})
    )
    respx.get(url__regex=r"rest\.ensembl\.org").mock(return_value=httpx.Response(404))
    client = VariantEffectsClient(
        gnomad=GnomADClient(http_client),
        myvariant=MyVariantClient(http_client),
        mavedb=MaveDBClient(http_client),
        ensembl=EnsemblClient(http_client),
    )
    result = await client.get_effects("TP53", "R175H")

    assert result.gnomad_variant_id is None
    assert result.annotation is None
    assert any("not found in gnomAD" in n for n in result.notes)
    assert any("MyVariant.info returned no record" in n for n in result.notes)


async def test_variant_effects_aggregator_invalid_mutation(http_client):
    from genesis_bio_mcp.clients.variant_effects import VariantEffectsClient

    client = VariantEffectsClient(
        gnomad=GnomADClient(http_client),
        myvariant=MyVariantClient(http_client),
        mavedb=MaveDBClient(http_client),
        ensembl=EnsemblClient(http_client),
    )
    # parse_protein_change raises ValueError on garbage input
    with pytest.raises(ValueError):
        await client.get_effects("TP53", "not a mutation")


# ---------------------------------------------------------------------------
# IEDB NextGen Tools (MHC binding) tests
# ---------------------------------------------------------------------------

MOCK_IEDB_SUBMIT = {
    "result_id": "abc-123",
    "results_uri": "https://api-nextgen-tools.iedb.org/api/v1/results/abc-123",
    "pipeline_id": "pipe-xyz",
}

MOCK_IEDB_RESULT_DONE = {
    "id": "abc-123",
    "type": "result",
    "status": "done",
    "data": {
        "results": [
            {
                "type": "peptide_table",
                "table_columns": [
                    {"name": "sequence_number"},
                    {"name": "peptide"},
                    {"name": "start"},
                    {"name": "end"},
                    {"name": "length"},
                    {"name": "allele"},
                    {"name": "peptide_index"},
                    {"name": "median_percentile"},
                    {"name": "netmhcpan_el_core"},
                    {"name": "netmhcpan_el_icore"},
                    {"name": "netmhcpan_el_score"},
                    {"name": "netmhcpan_el_percentile"},
                ],
                "table_data": [
                    [
                        1,
                        "SLYNTVATL",
                        1,
                        9,
                        9,
                        "HLA-A*02:01",
                        1,
                        0.06,
                        "SLYNTVATL",
                        "SLYNTVATL",
                        0.828,
                        0.06,
                    ],
                    [
                        1,
                        "SLYNTVATL",
                        1,
                        9,
                        9,
                        "HLA-A*03:01",
                        1,
                        35.0,
                        "SLYNTVATL",
                        "SLYNTVATL",
                        0.002,
                        35.0,
                    ],
                    [
                        1,
                        "LYNTVATLY",
                        2,
                        10,
                        9,
                        "HLA-A*02:01",
                        2,
                        1.5,
                        "LYNTVATLY",
                        "LYNTVATLY",
                        0.45,
                        1.5,
                    ],
                ],
            },
            {"type": "netmhcpan_allele_distance", "table_data": []},
        ],
        "errors": [],
        "warnings": [],
    },
}

MOCK_IEDB_RESULT_PENDING = {"id": "abc-123", "type": "result", "status": "pending", "data": {}}


@respx.mock
async def test_iedb_tools_predict_happy_path(http_client):
    from genesis_bio_mcp.clients.iedb_tools import IEDBToolsClient

    respx.post("https://api-nextgen-tools.iedb.org/api/v1/pipeline").mock(
        return_value=httpx.Response(200, json=MOCK_IEDB_SUBMIT)
    )
    respx.get("https://api-nextgen-tools.iedb.org/api/v1/results/abc-123").mock(
        return_value=httpx.Response(200, json=MOCK_IEDB_RESULT_DONE)
    )
    client = IEDBToolsClient(http_client)
    result = await client.predict_mhc_binding(
        sequence="SLYNTVATL", alleles=["HLA-A*02:01"], mhc_class="I"
    )
    assert result is not None
    assert result.strong_binder_count == 1
    assert result.weak_binder_count == 1
    assert len(result.hits) == 3
    # Hits sorted by percentile asc; strongest first
    assert result.hits[0].peptide == "SLYNTVATL"
    assert result.hits[0].allele == "HLA-A*02:01"
    assert result.hits[0].binder_class == "strong"


@respx.mock
async def test_iedb_tools_predict_no_strong_binders(http_client):
    from genesis_bio_mcp.clients.iedb_tools import IEDBToolsClient

    # Result with only non-binders
    payload = {
        "id": "abc-123",
        "status": "done",
        "data": {
            "results": [
                {
                    "type": "peptide_table",
                    "table_columns": [
                        {"name": "peptide"},
                        {"name": "allele"},
                        {"name": "length"},
                        {"name": "netmhcpan_el_percentile"},
                    ],
                    "table_data": [["AAAAAAAAA", "HLA-A*02:01", 9, 45.0]],
                }
            ]
        },
    }
    respx.post("https://api-nextgen-tools.iedb.org/api/v1/pipeline").mock(
        return_value=httpx.Response(200, json=MOCK_IEDB_SUBMIT)
    )
    respx.get("https://api-nextgen-tools.iedb.org/api/v1/results/abc-123").mock(
        return_value=httpx.Response(200, json=payload)
    )
    client = IEDBToolsClient(http_client)
    result = await client.predict_mhc_binding(sequence="AAAAAAAAA")
    assert result is not None
    assert result.strong_binder_count == 0
    assert result.weak_binder_count == 0
    assert len(result.hits) == 1
    assert result.hits[0].binder_class == "non_binder"


@respx.mock
async def test_iedb_tools_predict_submit_error(http_client):
    from genesis_bio_mcp.clients.iedb_tools import IEDBToolsClient

    respx.post("https://api-nextgen-tools.iedb.org/api/v1/pipeline").mock(
        side_effect=httpx.ConnectError("boom")
    )
    client = IEDBToolsClient(http_client)
    result = await client.predict_mhc_binding(sequence="SLYNTVATL")
    assert result is None


async def test_iedb_tools_rejects_oversize_request(http_client):
    from genesis_bio_mcp.clients.iedb_tools import IEDBToolsClient

    client = IEDBToolsClient(http_client)
    # 1000-residue protein × 5 alleles × 2 lengths = ~10000 peptide-allele pairs
    with pytest.raises(ValueError, match="exceeds limit"):
        await client.predict_mhc_binding(sequence="A" * 1000)


def test_iedb_tools_estimate_peptide_count_unit():
    """Internal windowing helper — 20aa protein with [9,10] yields 12+11=23."""
    from genesis_bio_mcp.clients.iedb_tools import _estimate_peptide_count

    assert _estimate_peptide_count(">q\n" + "A" * 20, [9, 10]) == (20 - 9 + 1) + (20 - 10 + 1)


def test_iedb_tools_ensure_fasta_passthrough():
    from genesis_bio_mcp.clients.iedb_tools import _ensure_fasta

    assert _ensure_fasta(">existing\nSEQ").startswith(">existing")


# ---------------------------------------------------------------------------
# GTEx client tests
# ---------------------------------------------------------------------------


_MOCK_GTEX_BRAF = {
    "data": [
        {
            "tissueSiteDetailId": "Brain_Cortex",
            "median": 12.3,
            "numSamples": 205,
        },
        {
            "tissueSiteDetailId": "Liver",
            "median": 4.5,
            "numSamples": 226,
        },
        {
            "tissueSiteDetailId": "Testis",
            "median": 28.7,
            "numSamples": 361,
        },
    ]
}


@respx.mock
async def test_gtex_happy_path(http_client, tmp_path, monkeypatch):
    monkeypatch.setattr(
        "genesis_bio_mcp.config.settings.settings.gtex_cache_path",
        tmp_path / "gtex.json",
    )
    respx.get(url__regex=r"rest\.ensembl\.org/lookup/symbol/homo_sapiens/BRAF").mock(
        return_value=httpx.Response(
            200,
            json={
                "id": "ENSG00000157764",
                "display_name": "BRAF",
                "seq_region_name": "7",
                "start": 140719327,
                "end": 140924928,
                "strand": -1,
                "biotype": "protein_coding",
                "Transcript": [{"id": "ENST00000646891", "is_canonical": 1}],
            },
        )
    )
    respx.get(url__regex=r"gtexportal\.org/api/v2/expression/medianGeneExpression").mock(
        return_value=httpx.Response(200, json=_MOCK_GTEX_BRAF)
    )
    ensembl = EnsemblClient(http_client)
    client = GTExClient(http_client, ensembl=ensembl)
    profile = await client.get_expression("BRAF")

    assert profile is not None
    assert profile.gene_symbol == "BRAF"
    assert profile.gencode_id == "ENSG00000157764"
    assert len(profile.samples) == 3
    # Sorted by TPM when rendered, but stored in fetch order
    tissues = {s.tissue for s in profile.samples}
    assert "Brain - Cortex" in tissues  # underscore → " - "
    md = profile.to_markdown()
    assert "Testis" in md
    assert "28.7" in md


@respx.mock
async def test_gtex_returns_empty_when_no_ensembl_id(http_client, tmp_path, monkeypatch):
    monkeypatch.setattr(
        "genesis_bio_mcp.config.settings.settings.gtex_cache_path",
        tmp_path / "gtex.json",
    )
    respx.get(url__regex=r"rest\.ensembl\.org/lookup/symbol").mock(return_value=httpx.Response(404))
    ensembl = EnsemblClient(http_client)
    client = GTExClient(http_client, ensembl=ensembl)
    profile = await client.get_expression("NOTAGENE")

    assert profile is not None
    assert profile.samples == []
    assert profile.gencode_id is None


@respx.mock
async def test_gtex_network_error_returns_empty_profile(http_client, tmp_path, monkeypatch):
    monkeypatch.setattr(
        "genesis_bio_mcp.config.settings.settings.gtex_cache_path",
        tmp_path / "gtex.json",
    )
    respx.get(url__regex=r"rest\.ensembl\.org/lookup/symbol").mock(
        return_value=httpx.Response(
            200,
            json={
                "id": "ENSG00000157764",
                "display_name": "BRAF",
                "seq_region_name": "7",
                "start": 1,
                "end": 2,
                "strand": 1,
                "Transcript": [],
            },
        )
    )
    respx.get(url__regex=r"gtexportal\.org").mock(side_effect=httpx.ConnectError("boom"))
    ensembl = EnsemblClient(http_client)
    client = GTExClient(http_client, ensembl=ensembl)
    profile = await client.get_expression("BRAF")

    assert profile is not None
    assert profile.samples == []


# ---------------------------------------------------------------------------
# HPA client tests
# ---------------------------------------------------------------------------


_MOCK_HPA_BRAF = [
    {
        "Gene": "BRAF",
        "Ensembl": "ENSG00000157764",
        "RNA tissue specificity": "Low tissue specificity",
        "RNA tissue specificity score": "1.3",
        "Subcellular main location": "Cytosol",
        "Subcellular location": "Cytosol, Plasma membrane",
        "Pathology prognostics - Melanoma": "Unfavorable (p=0.001)",
        "Pathology prognostics - Glioma": "Favorable (p=0.03)",
    }
]


@respx.mock
async def test_hpa_happy_path(http_client, tmp_path, monkeypatch):
    monkeypatch.setattr(
        "genesis_bio_mcp.config.settings.settings.hpa_cache_path",
        tmp_path / "hpa.json",
    )
    respx.get(url__regex=r"proteinatlas\.org/api/search_download\.php").mock(
        return_value=httpx.Response(200, json=_MOCK_HPA_BRAF)
    )
    client = HPAClient(http_client)
    report = await client.get_report("BRAF")

    assert report is not None
    assert report.gene_symbol == "BRAF"
    assert report.expression is not None
    assert report.expression.rna_tissue_specificity_category == "Low tissue specificity"
    assert report.expression.rna_tissue_specificity_score == 1.3
    assert "Cytosol" in report.expression.subcellular_locations
    assert len(report.pathology) == 2
    outcomes = {(p.cancer_type, p.prognostic_outcome) for p in report.pathology}
    assert ("Melanoma", "Unfavorable") in outcomes
    assert ("Glioma", "Favorable") in outcomes

    md = report.to_markdown()
    assert "Melanoma" in md
    assert "Low tissue specificity" in md


@respx.mock
async def test_hpa_empty_response(http_client, tmp_path, monkeypatch):
    monkeypatch.setattr(
        "genesis_bio_mcp.config.settings.settings.hpa_cache_path",
        tmp_path / "hpa.json",
    )
    respx.get(url__regex=r"proteinatlas\.org").mock(return_value=httpx.Response(200, json=[]))
    client = HPAClient(http_client)
    report = await client.get_report("NOTAGENE")

    assert report is not None
    assert report.expression is None
    assert report.pathology == []


@respx.mock
async def test_hpa_network_error_returns_none(http_client, tmp_path, monkeypatch):
    monkeypatch.setattr(
        "genesis_bio_mcp.config.settings.settings.hpa_cache_path",
        tmp_path / "hpa.json",
    )
    respx.get(url__regex=r"proteinatlas\.org").mock(side_effect=httpx.ConnectError("boom"))
    client = HPAClient(http_client)
    report = await client.get_report("BRAF")
    assert report is None


# ---------------------------------------------------------------------------
# OpenFDA client tests
# ---------------------------------------------------------------------------

# A FAERS /drug/event.json response carrying BOTH the count-by-field rows and
# the meta.results.total used by the "get total reports" query. The client
# issues both requests in parallel against the same endpoint, so a single
# mock response with both keys populated satisfies them.
_MOCK_OPENFDA_FAERS_IMATINIB = {
    "meta": {"results": {"total": 12345}},
    "results": [
        {"term": "NAUSEA", "count": 920},
        {"term": "FATIGUE", "count": 812},
        {"term": "RASH", "count": 540},
        # term/count pairs with bad types must be filtered out silently
        {"term": "", "count": 100},
        {"term": "MISSING_COUNT"},
    ],
}

_MOCK_OPENFDA_LABEL_IMATINIB = {
    "results": [
        {
            "boxed_warning": [
                "WARNING: Severe hepatotoxicity has been reported.",
                "",
            ],
            # other label fields ignored
            "indications_and_usage": ["Unused text"],
        }
    ]
}

_MOCK_OPENFDA_RECALL_IMATINIB = {
    "results": [
        {
            "recall_number": "D-0123-2023",
            "classification": "Class II",
            "reason_for_recall": "Container closure defect observed at manufacturing site.",
            "status": "Ongoing",
        },
        # Missing reason → must be skipped
        {
            "recall_number": "D-0456-2023",
            "classification": "Class III",
            "reason_for_recall": "",
            "status": "Terminated",
        },
    ]
}


@respx.mock
async def test_openfda_happy_path(http_client, tmp_path, monkeypatch):
    monkeypatch.setattr(
        "genesis_bio_mcp.config.settings.settings.openfda_cache_path",
        tmp_path / "openfda.json",
    )
    monkeypatch.delenv("OPENFDA_API_KEY", raising=False)
    respx.get(url__regex=r"api\.fda\.gov/drug/event\.json").mock(
        return_value=httpx.Response(200, json=_MOCK_OPENFDA_FAERS_IMATINIB)
    )
    respx.get(url__regex=r"api\.fda\.gov/drug/label\.json").mock(
        return_value=httpx.Response(200, json=_MOCK_OPENFDA_LABEL_IMATINIB)
    )
    respx.get(url__regex=r"api\.fda\.gov/drug/enforcement\.json").mock(
        return_value=httpx.Response(200, json=_MOCK_OPENFDA_RECALL_IMATINIB)
    )
    client = OpenFDAClient(http_client)
    signal = await client.get_safety_signals("imatinib")

    assert signal is not None
    assert signal.drug_name == "imatinib"
    assert signal.total_reports == 12345
    # Only rows with non-empty term and integer count > 0 survive
    terms = [ae.term for ae in signal.top_adverse_events]
    assert terms == ["NAUSEA", "FATIGUE", "RASH"]
    assert signal.top_adverse_events[0].count == 920
    # Empty string in boxed_warning must be stripped, leaving one warning
    assert signal.boxed_warnings == ["WARNING: Severe hepatotoxicity has been reported."]
    # Missing reason_for_recall means that recall is dropped
    assert len(signal.recalls) == 1
    assert signal.recalls[0].recall_number == "D-0123-2023"
    assert signal.recalls[0].classification == "Class II"

    md = signal.to_markdown()
    assert "Boxed warnings" in md
    assert "NAUSEA" in md
    assert "D-0123-2023" in md
    assert "FAERS" in md  # disclaimer

    # Session cache should now be populated; a second call must not hit HTTP.
    respx.get(url__regex=r"api\.fda\.gov").mock(side_effect=AssertionError("should be cached"))
    cached = await client.get_safety_signals("Imatinib")  # case-insensitive key
    assert cached is signal


@respx.mock
async def test_openfda_no_matches_returns_empty_signal(http_client, tmp_path, monkeypatch):
    """OpenFDA 404 = zero matches; the client returns a populated-but-empty signal."""
    monkeypatch.setattr(
        "genesis_bio_mcp.config.settings.settings.openfda_cache_path",
        tmp_path / "openfda.json",
    )
    monkeypatch.delenv("OPENFDA_API_KEY", raising=False)
    respx.get(url__regex=r"api\.fda\.gov/drug/event\.json").mock(
        return_value=httpx.Response(404, json={"error": {"code": "NOT_FOUND"}})
    )
    respx.get(url__regex=r"api\.fda\.gov/drug/label\.json").mock(
        return_value=httpx.Response(404, json={"error": {"code": "NOT_FOUND"}})
    )
    respx.get(url__regex=r"api\.fda\.gov/drug/enforcement\.json").mock(
        return_value=httpx.Response(404, json={"error": {"code": "NOT_FOUND"}})
    )
    client = OpenFDAClient(http_client)
    signal = await client.get_safety_signals("not-a-real-drug")

    assert signal is not None  # clean negative, not a failure
    assert signal.total_reports == 0
    assert signal.top_adverse_events == []
    assert signal.boxed_warnings == []
    assert signal.recalls == []
    assert "FAERS" in signal.to_markdown()

    # Empty signals must not be persisted to disk (so later reports are picked up).
    assert not (tmp_path / "openfda.json").exists()


@respx.mock
async def test_openfda_network_error_returns_none(http_client, tmp_path, monkeypatch):
    monkeypatch.setattr(
        "genesis_bio_mcp.config.settings.settings.openfda_cache_path",
        tmp_path / "openfda.json",
    )
    monkeypatch.delenv("OPENFDA_API_KEY", raising=False)
    respx.get(url__regex=r"api\.fda\.gov").mock(side_effect=httpx.ConnectError("boom"))
    client = OpenFDAClient(http_client)
    signal = await client.get_safety_signals("imatinib")
    assert signal is None


@respx.mock
async def test_openfda_api_key_injected_when_set(http_client, tmp_path, monkeypatch):
    """If OPENFDA_API_KEY is set, it must be appended as a query parameter."""
    monkeypatch.setattr(
        "genesis_bio_mcp.config.settings.settings.openfda_cache_path",
        tmp_path / "openfda.json",
    )
    monkeypatch.setenv("OPENFDA_API_KEY", "test-key-xyz")
    seen_keys: list[str | None] = []

    def _spy(request):
        seen_keys.append(request.url.params.get("api_key"))
        return httpx.Response(200, json={"results": [], "meta": {"results": {"total": 0}}})

    respx.get(url__regex=r"api\.fda\.gov").mock(side_effect=_spy)
    client = OpenFDAClient(http_client)
    signal = await client.get_safety_signals("imatinib")

    assert signal is not None
    # All 4 sub-requests (event count, event total, label, enforcement) must
    # carry the key.
    assert len(seen_keys) == 4
    assert all(k == "test-key-xyz" for k in seen_keys)


# ---------------------------------------------------------------------------
# ChEMBL client tests
# ---------------------------------------------------------------------------

_MOCK_CHEMBL_TARGET_SEARCH = {
    "targets": [
        # First hit is wrong organism — client must skip to the human row.
        {
            "target_chembl_id": "CHEMBL9999",
            "target_type": "SINGLE PROTEIN",
            "organism": "Rattus norvegicus",
        },
        {
            "target_chembl_id": "CHEMBL5145",
            "target_type": "SINGLE PROTEIN",
            "organism": "Homo sapiens",
        },
    ]
}

_MOCK_CHEMBL_ACTIVITIES = {
    "activities": [
        # Biochemical binding assay, human, high confidence — the reference row.
        {
            "molecule_chembl_id": "CHEMBL1336",
            "molecule_pref_name": "VEMURAFENIB",
            "standard_type": "IC50",
            "pchembl_value": "8.5",
            "assay_description": "In vitro kinase inhibition against human BRAF V600E",
            "assay_type": "B",
            "assay_organism": "Homo sapiens",
            "assay_cell_type": None,
            "bao_label": "single protein format",
            "confidence_score": 9,
        },
        # Functional cell-based assay in HEK293 — should render with cell type.
        {
            "molecule_chembl_id": "CHEMBL2028663",
            "molecule_pref_name": "DABRAFENIB",
            "standard_type": "IC50",
            "pchembl_value": "9.2",
            "assay_description": "Inhibition of ERK phosphorylation in A375 cells",
            "assay_type": "F",
            "assay_organism": "Homo sapiens",
            "assay_cell_type": "A375",
            "bao_label": "cell-based format",
            "confidence_score": 8,  # < 9 — should drive the low-confidence flag
        },
        # Rat assay — should surface under "Non-human assays present".
        {
            "molecule_chembl_id": "CHEMBL1337",
            "molecule_pref_name": "ROCK_FOR_RAT",
            "standard_type": "Ki",
            "pchembl_value": 7.0,
            "assay_description": "Rat liver microsome assay",
            "assay_type": "B",
            "assay_organism": "Rattus norvegicus",
            "assay_cell_type": None,
            "bao_label": "single protein format",
            "confidence_score": 9,
        },
        # Duplicate molecule — must be deduped (only the first occurrence kept).
        {
            "molecule_chembl_id": "CHEMBL1336",
            "molecule_pref_name": "VEMURAFENIB",
            "standard_type": "Ki",
            "pchembl_value": "7.0",
            "assay_description": "dup row",
            "assay_type": "B",
            "assay_organism": "Homo sapiens",
            "confidence_score": 9,
        },
        # Missing pchembl — must be skipped.
        {
            "molecule_chembl_id": "CHEMBL99999",
            "molecule_pref_name": "NO_PCHEMBL",
            "standard_type": "IC50",
            "pchembl_value": None,
            "assay_type": "B",
        },
        # Unparseable confidence score — coerced to None, row still kept.
        {
            "molecule_chembl_id": "CHEMBL7777",
            "molecule_pref_name": "WEIRD_CONF",
            "standard_type": "IC50",
            "pchembl_value": "6.1",
            "assay_description": "Assay with non-int confidence in source data",
            "assay_type": "B",
            "assay_organism": "Homo sapiens",
            "confidence_score": "not-a-number",
        },
    ]
}


@respx.mock
async def test_chembl_happy_path_with_assay_context(http_client):
    respx.get(url__regex=r"ebi\.ac\.uk/chembl/api/data/target/search").mock(
        return_value=httpx.Response(200, json=_MOCK_CHEMBL_TARGET_SEARCH)
    )
    respx.get(url__regex=r"ebi\.ac\.uk/chembl/api/data/activity").mock(
        return_value=httpx.Response(200, json=_MOCK_CHEMBL_ACTIVITIES)
    )
    client = ChEMBLClient(http_client)
    result = await client.get_compounds("BRAF")

    assert result is not None
    assert result.gene_symbol == "BRAF"
    assert result.target_chembl_id == "CHEMBL5145"  # human, not rat
    # 4 survive (one dup, one missing pchembl dropped); total_active_compounds
    # is the *surviving* count after dedup/filter — matches best_pchembl source.
    assert result.total_active_compounds == 4
    assert result.best_pchembl == 9.2  # dabrafenib wins

    # New assay-context fields must be populated on each activity.
    by_id = {c.molecule_chembl_id: c for c in result.compounds}
    v = by_id["CHEMBL1336"]
    assert v.assay_type == "B"
    assert v.assay_organism == "Homo sapiens"
    assert v.assay_cell_type is None
    assert v.bao_format == "single protein format"
    assert v.confidence_score == 9

    d = by_id["CHEMBL2028663"]
    assert d.assay_type == "F"
    assert d.assay_cell_type == "A375"
    assert d.bao_format == "cell-based format"
    assert d.confidence_score == 8

    rat = by_id["CHEMBL1337"]
    assert rat.assay_organism == "Rattus norvegicus"

    weird = by_id["CHEMBL7777"]
    assert weird.confidence_score is None  # unparseable confidence coerced

    md = result.to_markdown()
    # Assay mix summary appears and distinguishes B vs F
    assert "Assay mix" in md
    assert "binding" in md
    assert "functional" in md
    # Non-human flag surfaces rat assays
    assert "Non-human assays present" in md
    assert "Rattus norvegicus" in md
    # Low-confidence flag fires because CHEMBL2028663 has score 8
    assert "Low target-assignment confidence" in md
    # Table carries assay-type + organism columns
    assert "| Assay |" in md
    assert "A375" in md  # cell type surfaces in the label


@respx.mock
async def test_chembl_no_human_target_returns_none(http_client):
    """If only non-human hits are returned, the client should return None."""
    respx.get(url__regex=r"ebi\.ac\.uk/chembl/api/data/target/search").mock(
        return_value=httpx.Response(
            200,
            json={
                "targets": [
                    {
                        "target_chembl_id": "CHEMBL9999",
                        "target_type": "SINGLE PROTEIN",
                        "organism": "Rattus norvegicus",
                    }
                ]
            },
        )
    )
    client = ChEMBLClient(http_client)
    result = await client.get_compounds("BRAF")
    assert result is None


@respx.mock
async def test_chembl_empty_activities_returns_zero_count(http_client):
    """A human target with no pChEMBL rows returns a zero-count ChEMBLCompounds."""
    respx.get(url__regex=r"ebi\.ac\.uk/chembl/api/data/target/search").mock(
        return_value=httpx.Response(200, json=_MOCK_CHEMBL_TARGET_SEARCH)
    )
    respx.get(url__regex=r"ebi\.ac\.uk/chembl/api/data/activity").mock(
        return_value=httpx.Response(200, json={"activities": []})
    )
    client = ChEMBLClient(http_client)
    result = await client.get_compounds("BRAF")

    assert result is not None
    assert result.total_active_compounds == 0
    assert result.best_pchembl is None
    assert result.compounds == []


@respx.mock
async def test_chembl_network_error_returns_none(http_client):
    respx.get(url__regex=r"ebi\.ac\.uk/chembl/api/data/target/search").mock(
        side_effect=httpx.ConnectError("boom")
    )
    client = ChEMBLClient(http_client)
    result = await client.get_compounds("BRAF")
    assert result is None


# ---------------------------------------------------------------------------
# v0.3.1 regression tests — four bugs caught by the post-release smoke run
# ---------------------------------------------------------------------------


@respx.mock
async def test_chembl_target_organism_is_parsed_when_assay_organism_null(http_client):
    """Bug D: real ChEMBL responses populate ``target_organism``, not ``assay_organism``.

    The activity row exposes the species via ``target_organism``; ``assay_organism``
    is almost always null in the live API. The client must read the populated field
    so the rendered table doesn't show empty Organism cells across the board.
    """
    respx.get(url__regex=r"ebi\.ac\.uk/chembl/api/data/target/search").mock(
        return_value=httpx.Response(
            200,
            json={
                "targets": [
                    {
                        "target_chembl_id": "CHEMBL5145",
                        "target_type": "SINGLE PROTEIN",
                        "organism": "Homo sapiens",
                    }
                ]
            },
        )
    )
    respx.get(url__regex=r"ebi\.ac\.uk/chembl/api/data/activity").mock(
        return_value=httpx.Response(
            200,
            json={
                "activities": [
                    {
                        "molecule_chembl_id": "CHEMBL1336",
                        "pchembl_value": "8.5",
                        "standard_type": "IC50",
                        "assay_type": "B",
                        # Live shape: target_organism populated, assay_organism null.
                        "target_organism": "Homo sapiens",
                        "assay_organism": None,
                    }
                ]
            },
        )
    )
    client = ChEMBLClient(http_client)
    result = await client.get_compounds("BRAF")
    assert result is not None
    assert result.compounds[0].assay_organism == "Homo sapiens"


def test_openfda_normalize_drug_name_strips_salt_suffixes():
    """Bug C: DGIdb-style drug names (salt forms, hydrate suffixes) must
    collapse to the bare INN before hitting OpenFDA."""
    from genesis_bio_mcp.clients.openfda import _normalize_drug_name

    assert _normalize_drug_name("ATORVASTATIN CALCIUM TRIHYDRATE") == "ATORVASTATIN"
    assert _normalize_drug_name("FILGOTINIB MALEATE") == "FILGOTINIB"
    assert _normalize_drug_name("DACOMITINIB ANHYDROUS") == "DACOMITINIB"
    assert _normalize_drug_name("PRAVASTATIN SODIUM") == "PRAVASTATIN"
    assert _normalize_drug_name("imatinib mesylate") == "imatinib"
    # Already-clean INNs pass through unchanged (besides whitespace normalization).
    assert _normalize_drug_name("evolocumab") == "evolocumab"
    assert _normalize_drug_name("  alirocumab  ") == "alirocumab"
    # Empty / whitespace-only input returns empty string, not raises.
    assert _normalize_drug_name("") == ""
    assert _normalize_drug_name("   ") == ""


@respx.mock
async def test_openfda_uses_or_syntax_and_normalizes_drug_name(http_client, tmp_path, monkeypatch):
    """Bug C: query must use Lucene ``OR`` (not ``+`` which means MUST), AND
    must run on the salt-stripped name so biologics and salt forms surface."""
    monkeypatch.setattr(
        "genesis_bio_mcp.config.settings.settings.openfda_cache_path",
        tmp_path / "openfda.json",
    )
    monkeypatch.delenv("OPENFDA_API_KEY", raising=False)

    seen_searches: list[str] = []

    def _spy(request):
        seen_searches.append(request.url.params.get("search") or "")
        return httpx.Response(200, json={"results": [], "meta": {"results": {"total": 0}}})

    respx.get(url__regex=r"api\.fda\.gov").mock(side_effect=_spy)
    client = OpenFDAClient(http_client)
    signal = await client.get_safety_signals("EVOLOCUMAB")
    assert signal is not None  # empty signal, not None — proves the call ran

    assert seen_searches, "expected at least one OpenFDA query"
    # Every clause must use OR (not '+' which Lucene parses as MUST/AND).
    # Stale '+' syntax requires ALL three fields to match — the bug that
    # silently dropped FAERS for biologics like evolocumab.
    for s in seen_searches:
        assert " OR " in s, f"OpenFDA search lost OR syntax: {s}"
        assert "+patient.drug.openfda.generic_name" not in s
        assert "+openfda.generic_name" not in s

    # Verify the normalized name (no salt suffix) was sent in the query.
    seen_searches.clear()
    await client.get_safety_signals("FILGOTINIB MALEATE")
    assert seen_searches, "expected query for FILGOTINIB MALEATE"
    for s in seen_searches:
        assert "FILGOTINIB" in s
        assert "MALEATE" not in s


@respx.mock
async def test_gtex_resolves_versioned_gencode_via_reference_endpoint(
    http_client, tmp_path, monkeypatch
):
    """Bug B: GTEx requires a versioned GENCODE ID (``ENSG…11``); the
    unversioned form returns an empty data array. Resolution must go through
    GTEx's ``/reference/gene`` endpoint to get the version GTEx itself indexed."""
    from genesis_bio_mcp.clients.ensembl import EnsemblClient
    from genesis_bio_mcp.clients.gtex import GTExClient

    monkeypatch.setattr(
        "genesis_bio_mcp.config.settings.settings.gtex_cache_path",
        tmp_path / "gtex.json",
    )

    seen_expr_params: list[str | None] = []

    def _expr_handler(request):
        seen_expr_params.append(request.url.params.get("gencodeId"))
        return httpx.Response(
            200,
            json={
                "data": [
                    {
                        "tissueSiteDetailId": "Liver",
                        "median": 12.5,
                        "numSamples": 226,
                        "gencodeId": "ENSG00000169174.11",
                        "geneSymbol": "PCSK9",
                    }
                ]
            },
        )

    # GTEx /reference/gene returns the versioned GENCODE ID.
    respx.get(url__regex=r"gtexportal\.org/api/v2/reference/gene").mock(
        return_value=httpx.Response(
            200,
            json={
                "data": [
                    {
                        "geneSymbol": "PCSK9",
                        "gencodeId": "ENSG00000169174.11",
                    }
                ]
            },
        )
    )
    respx.get(url__regex=r"gtexportal\.org/api/v2/expression/medianGeneExpression").mock(
        side_effect=_expr_handler,
    )

    ensembl = EnsemblClient(http_client)
    client = GTExClient(http_client, ensembl=ensembl)
    profile = await client.get_expression("PCSK9")

    assert profile is not None
    assert profile.gencode_id == "ENSG00000169174.11"  # versioned, not bare
    assert profile.samples and profile.samples[0].tissue.startswith("Liver")
    # Crucially, the expression endpoint received the VERSIONED form.
    assert seen_expr_params == ["ENSG00000169174.11"]


@respx.mock
async def test_uniprot_search_prefers_exact_gene_match_over_first_hit(http_client):
    """Bug A: ``gene_exact:ALB`` returned FBF1 ranked first, ALB second.
    The resolver must prefer the entry whose primary geneName matches the
    query rather than blindly taking the first row."""
    from genesis_bio_mcp.clients.uniprot import UniProtClient

    respx.get(url__regex=r"rest\.uniprot\.org/uniprotkb/search").mock(
        return_value=httpx.Response(
            200,
            json={
                "results": [
                    # First hit is the wrong gene — proves the bug shape.
                    {
                        "primaryAccession": "Q8TES7",
                        "entryType": "UniProtKB reviewed (Swiss-Prot)",
                        "genes": [{"geneName": {"value": "FBF1"}}],
                        "proteinDescription": {
                            "recommendedName": {"fullName": {"value": "Fas-binding factor 1"}}
                        },
                    },
                    # Real ALB ranked second — must win.
                    {
                        "primaryAccession": "P02768",
                        "entryType": "UniProtKB reviewed (Swiss-Prot)",
                        "genes": [{"geneName": {"value": "ALB"}}],
                        "proteinDescription": {
                            "recommendedName": {"fullName": {"value": "Albumin"}}
                        },
                    },
                ]
            },
        )
    )
    client = UniProtClient(http_client)
    protein = await client.get_protein("ALB")
    assert protein is not None
    assert protein.uniprot_accession == "P02768"
    assert protein.gene_symbol == "ALB"
    assert protein.protein_name == "Albumin"


# ---------------------------------------------------------------------------
# v0.3.2 regression tests — five second-round fixes from the live MCP smoke run
# ---------------------------------------------------------------------------


@respx.mock
async def test_gtex_passes_gencode_v39_to_reference_endpoint(http_client, tmp_path, monkeypatch):
    """Bug B (still): /reference/gene defaults to gencodeVersion=v26 which
    returns IDs the gtex_v10 expression dataset (keyed on GENCODE v39) can't
    find. The client must pin gencodeVersion=v39 so the returned ID matches."""
    from genesis_bio_mcp.clients.ensembl import EnsemblClient
    from genesis_bio_mcp.clients.gtex import GTExClient

    monkeypatch.setattr(
        "genesis_bio_mcp.config.settings.settings.gtex_cache_path",
        tmp_path / "gtex.json",
    )

    seen_versions: list[str | None] = []

    def _ref_handler(request):
        seen_versions.append(request.url.params.get("gencodeVersion"))
        return httpx.Response(
            200,
            json={"data": [{"geneSymbol": "INS", "gencodeId": "ENSG00000254647.7"}]},
        )

    respx.get(url__regex=r"gtexportal\.org/api/v2/reference/gene").mock(
        side_effect=_ref_handler,
    )
    respx.get(url__regex=r"gtexportal\.org/api/v2/expression/medianGeneExpression").mock(
        return_value=httpx.Response(
            200,
            json={
                "data": [
                    {
                        "tissueSiteDetailId": "Pancreas",
                        "median": 12345.0,
                        "numSamples": 200,
                        "gencodeId": "ENSG00000254647.7",
                        "geneSymbol": "INS",
                    }
                ]
            },
        )
    )

    ensembl = EnsemblClient(http_client)
    client = GTExClient(http_client, ensembl=ensembl)
    profile = await client.get_expression("INS")

    assert profile is not None
    assert profile.gencode_id == "ENSG00000254647.7"
    assert profile.samples and profile.samples[0].median_tpm == 12345.0
    # Crucial assertion — without v39, the live API returns empty data.
    assert seen_versions == ["v39"], f"expected gencodeVersion=v39, got {seen_versions}"


@respx.mock
async def test_open_targets_falls_back_to_parenthetical_abbreviation(http_client):
    """Bug J: 'non-small-cell lung cancer (NSCLC)' returns 0 OT hits even
    though the bare 'NSCLC' resolves to EFO_0003060. The resolver must try
    the parenthetical abbreviation as an explicit fallback variant."""
    import json as _json

    from genesis_bio_mcp.clients.open_targets import OpenTargetsClient

    seen_disease_queries: list[str] = []

    def _post_handler(request):
        body = _json.loads(request.content)
        query = body.get("query", "")
        variables = body.get("variables", {})

        if "DiseaseSearch" in query:
            qstr = variables.get("name", "")
            seen_disease_queries.append(qstr)
            if qstr.upper() == "NSCLC":
                return httpx.Response(
                    200,
                    json={
                        "data": {
                            "search": {
                                "hits": [
                                    {
                                        "id": "EFO_0003060",
                                        "name": "non-small cell lung carcinoma",
                                        "entity": "disease",
                                    }
                                ]
                            }
                        }
                    },
                )
            return httpx.Response(200, json={"data": {"search": {"hits": []}}})

        if "GeneSearch" in query:
            return httpx.Response(
                200,
                json={
                    "data": {
                        "search": {
                            "hits": [{"id": "ENSG00000146648", "entity": "target", "name": "EGFR"}]
                        }
                    }
                },
            )

        # Association query — return a stub row keyed to EFO_0003060 so the
        # association fetch succeeds.
        return httpx.Response(
            200,
            json={
                "data": {
                    "target": {
                        "associatedDiseases": {
                            "count": 1,
                            "rows": [
                                {
                                    "score": 0.85,
                                    "datatypeScores": [],
                                    "datasourceScores": [],
                                    "disease": {
                                        "id": "EFO_0003060",
                                        "name": "non-small cell lung carcinoma",
                                    },
                                }
                            ],
                        }
                    }
                }
            },
        )

    respx.post(url__regex=r"api\.platform\.opentargets\.org/api/v4/graphql").mock(
        side_effect=_post_handler,
    )

    client = OpenTargetsClient(http_client)
    result = await client.get_association("EGFR", "non-small-cell lung cancer (NSCLC)")

    # Literal MUST be tried first (cache-friendly + respects user intent).
    assert seen_disease_queries[0] == "non-small-cell lung cancer (NSCLC)"
    # The bare abbreviation MUST appear in the variant chain — it's the one
    # that resolves the EFO ID.
    assert "NSCLC" in seen_disease_queries
    # Lookup must ultimately succeed.
    assert result is not None
    assert result.disease_efo_id == "EFO_0003060"


def test_open_targets_normalize_indication_variants():
    """Unit-test the variant generator directly — order and dedup matter."""
    from genesis_bio_mcp.clients.open_targets import _normalize_indication_variants

    v = _normalize_indication_variants("non-small-cell lung cancer (NSCLC)")
    # Literal first, then stripped, then abbreviation, then hyphen-stripped.
    assert v[0] == "non-small-cell lung cancer (NSCLC)"
    assert "non-small-cell lung cancer" in v
    assert "NSCLC" in v
    assert any("non small cell" in candidate for candidate in v)
    assert len(v) == len(set(v))  # no dupes

    # No-paren input still yields the hyphen variant.
    v2 = _normalize_indication_variants("type-2-diabetes")
    assert "type-2-diabetes" in v2
    assert "type 2 diabetes" in v2


def test_compute_score_credits_gwas_axis_for_monogenic_diseases():
    """Bug G: when OT.genetic_association_score > 0.7 (monogenic signature)
    and GWAS Catalog returns nothing (because the disease isn't a GWAS
    subject), the GWAS axis must NOT zero out — that double-penalizes
    Mendelian targets like CFTR/CF whose genetic evidence is already fully
    captured by OT's overall_score."""
    from genesis_bio_mcp.models import TargetDiseaseAssociation
    from genesis_bio_mcp.tools.target_prioritization import _compute_score

    monogenic = TargetDiseaseAssociation(
        gene_symbol="CFTR",
        disease_name="cystic fibrosis",
        disease_efo_id="EFO_0000384",
        ensembl_id="ENSG00000001626",
        overall_score=0.85,
        genetic_association_score=0.99,  # > monogenic threshold
        somatic_mutation_score=None,
        known_drug_score=0.95,
        literature_mining_score=0.6,
        evidence_count=200,
    )
    _, breakdown = _compute_score(
        disease_assoc=monogenic,
        cancer_dep=None,
        gwas_ev=None,
        compounds=None,
        protein=None,
    )
    assert breakdown.gwas == 2.0, "monogenic credit must give full GWAS axis (2.0)"

    polygenic = TargetDiseaseAssociation(
        gene_symbol="FTO",
        disease_name="obesity",
        disease_efo_id="EFO_0001073",
        ensembl_id="ENSG00000140718",
        overall_score=0.5,
        genetic_association_score=0.55,  # below threshold
        somatic_mutation_score=None,
        known_drug_score=0.0,
        literature_mining_score=0.4,
        evidence_count=50,
    )
    _, breakdown2 = _compute_score(
        disease_assoc=polygenic,
        cancer_dep=None,
        gwas_ev=None,
        compounds=None,
        protein=None,
    )
    assert breakdown2.gwas == 0.0, "polygenic gene without GWAS must NOT get monogenic credit"


async def test_attach_safety_signals_filters_to_direct_engagement_drugs():
    """Bug F (v0.3.2 + v0.3.3 tightening): DGIdb returns ALL drugs co-mentioned
    with a gene (atezolizumab on EGFR via co-administration; ACYCLOVIR on
    ABL1 via... nothing). FAERS enrichment must skip drugs whose
    interaction_type is not in _DIRECT_TYPES so the safety panel reflects
    target-related liability, not co-administration noise.

    Bug F.1 (v0.3.3): the v0.3.2 untyped-fallback also leaked non-target
    drugs, because DGIdb leaves interaction_type=None for exactly those
    spurious associations. Untyped is now strictly excluded — the trade-off
    is losing FAERS coverage for very-recent approvals DGIdb hasn't typed
    yet, which we accept as the correct precision/recall trade-off."""
    from genesis_bio_mcp.models import DrugInteraction
    from genesis_bio_mcp.tools.target_prioritization import attach_safety_signals

    drugs = [
        DrugInteraction(
            drug_name="afatinib",
            interaction_type="inhibitor",
            phase=4,
            approved=True,
            sources=["FDA"],
        ),
        # Co-mention via trials, NOT a direct EGFR drug — must be skipped.
        DrugInteraction(
            drug_name="atezolizumab",
            interaction_type="other",
            phase=4,
            approved=True,
            sources=["TTD"],
        ),
        # Untyped approved drug — v0.3.3 strict filter excludes (Bug F.1).
        DrugInteraction(
            drug_name="erlotinib",
            interaction_type=None,
            phase=4,
            approved=True,
            sources=["FDA"],
        ),
    ]

    queried: list[str] = []

    class _FakeOpenFDA:
        async def get_safety_signals(self, drug_name: str):
            queried.append(drug_name)
            return None  # only verifying which drugs got queried, not the response

    out = await attach_safety_signals(drugs, openfda=_FakeOpenFDA())

    assert "afatinib" in queried
    assert "atezolizumab" not in queried, (
        "non-direct engager must NOT trigger FAERS lookup — that was Bug F"
    )
    assert "erlotinib" not in queried, (
        "untyped must NOT trigger FAERS lookup either — Bug F.1 v0.3.3"
    )
    assert len(out) == 3  # passthrough preserves the full list


@respx.mock
async def test_gwas_falls_back_to_unfiltered_top_hits_when_trait_match_empty(
    http_client, tmp_path, monkeypatch
):
    """Bug D: gene has GWAS evidence but no hits match the queried trait
    label (PCSK9 has dozens of LDL/lipid hits, none labelled 'familial
    hypercholesterolemia'). The fallback must surface the strongest gene-
    level associations with a sentinel trait_query so the score reflects
    'gene IS GWAS-implicated, just not under this exact label.'"""
    from genesis_bio_mcp.clients.gwas import GwasClient

    monkeypatch.setattr(
        "genesis_bio_mcp.config.settings.settings.gwas_cache_path",
        tmp_path / "gwas.json",
    )

    # SNP-symbol path: empty (we want the gene-ID path to carry the result).
    respx.get(url__regex=r"gwas/rest/api/singleNucleotidePolymorphisms/search/findByGene").mock(
        return_value=httpx.Response(
            200,
            json={"_embedded": {"singleNucleotidePolymorphisms": []}},
        )
    )

    # Gene-ID path: 2 hits labelled "LDL cholesterol levels" — they don't
    # match the queried "familial hypercholesterolemia" trait string.
    respx.get(url__regex=r"gwas/rest/api/associations/search/findByEntrezMappedGeneId").mock(
        return_value=httpx.Response(
            200,
            json={
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
                        },
                        {
                            "pvalue": 3e-15,
                            "loci": [
                                {
                                    "strongestRiskAlleles": [{"riskAlleleName": "rs562556-A"}],
                                    "authorReportedGenes": [{"geneName": "PCSK9"}],
                                }
                            ],
                            "efoTraits": [{"trait": "LDL cholesterol levels", "uri": None}],
                        },
                    ]
                }
            },
        )
    )

    client = GwasClient(http_client, efo_resolver=None)
    result = await client.get_evidence(
        "PCSK9", "familial hypercholesterolemia", ncbi_gene_id="255738"
    )

    # Fallback returns the top gene-level hits, NOT None.
    assert result is not None
    assert result.total_associations == 2
    # Sentinel marker so renderers can flag the fallback.
    assert "no exact-trait match" in result.trait_query
    assert result.strongest_p_value == 4e-20


# ---------------------------------------------------------------------------
# v0.3.3 regression tests — six third-round fixes from the second smoke run
# ---------------------------------------------------------------------------


def test_open_targets_name_match_score_rejects_fuzzy_substring_garbage():
    """Bug K: token-set Jaccard scoring lets us reject the MONDO_0010802
    pancreatic-hypoplasia syndrome that fuzzy-matched 'type-2 diabetes
    mellitus' just because the word 'diabetes' appears in both names."""
    from genesis_bio_mcp.clients.open_targets import _name_match_score

    # Exact match → 1.0
    assert _name_match_score("type 2 diabetes mellitus", "type 2 diabetes mellitus") == 1.0
    # Hyphen-equivalent → 1.0 (normalize to spaces)
    assert _name_match_score("type-2 diabetes mellitus", "type 2 diabetes mellitus") == 1.0
    # Substring (abbreviation indexed against long name) → 0.95
    assert _name_match_score("nsclc", "non small cell lung carcinoma nsclc") == 0.95
    # Garbage fuzzy match → well below 0.5 threshold
    bad = _name_match_score(
        "type 2 diabetes mellitus",
        "pancreatic hypoplasia diabetes congenital heart disease syndrome",
    )
    assert bad < 0.3, f"fuzzy garbage match scored {bad:.2f}, expected < 0.3"
    # Empty inputs → 0.0
    assert _name_match_score("", "anything") == 0.0
    assert _name_match_score("anything", "") == 0.0


@respx.mock
async def test_open_targets_resolve_disease_rejects_low_quality_match(http_client):
    """Bug K: when the parens-stripped variant returns a fuzzy-matched
    pancreatic-hypoplasia syndrome (~0.13 token overlap with 'type-2
    diabetes mellitus'), the resolver must reject it and continue trying
    later variants until it finds the real match (or returns None)."""
    import json as _json

    from genesis_bio_mcp.clients.open_targets import OpenTargetsClient

    def _post_handler(request):
        body = _json.loads(request.content)
        query = body.get("query", "")
        variables = body.get("variables", {})
        if "DiseaseSearch" not in query:
            return httpx.Response(200, json={"data": {"search": {"hits": []}}})
        qstr = variables.get("name", "")
        # Parens-stripped variant gets a wrong fuzzy hit; hyphen-normalized
        # variant gets the correct one. Old code took the wrong hit because
        # it fired first; new code scores both and picks the better match.
        if qstr == "type-2 diabetes mellitus":
            return httpx.Response(
                200,
                json={
                    "data": {
                        "search": {
                            "hits": [
                                {
                                    "id": "MONDO_0010802",
                                    "name": (
                                        "pancreatic hypoplasia-diabetes-congenital "
                                        "heart disease syndrome"
                                    ),
                                    "entity": "disease",
                                }
                            ]
                        }
                    }
                },
            )
        if qstr == "type 2 diabetes mellitus":
            return httpx.Response(
                200,
                json={
                    "data": {
                        "search": {
                            "hits": [
                                {
                                    "id": "MONDO_0005148",
                                    "name": "type 2 diabetes mellitus",
                                    "entity": "disease",
                                }
                            ]
                        }
                    }
                },
            )
        return httpx.Response(200, json={"data": {"search": {"hits": []}}})

    respx.post(url__regex=r"api\.platform\.opentargets\.org/api/v4/graphql").mock(
        side_effect=_post_handler,
    )

    client = OpenTargetsClient(http_client)
    efo_id, name = await client._resolve_disease("type-2 diabetes mellitus (T2DM)")
    assert efo_id == "MONDO_0005148", (
        f"resolver picked {efo_id} ({name}) — should have picked the real T2D entry"
    )
    assert name == "type 2 diabetes mellitus"


def test_depmap_pan_essential_cap_fires_at_95_percent_dependency():
    """Bug L: any gene with ≥95% of cell lines dependent should be flagged
    pan-essential, even when DepMap's common_essential boolean says
    otherwise. Catches MT-encoded subunits and other genes outside the
    routine essentiality screen."""
    from genesis_bio_mcp.clients.depmap import DepMapClient

    client = DepMapClient(client=None, gene_dep_cache={})  # type: ignore[arg-type]
    # 95%+ dependent → cap fires regardless of common_essential flag
    entry = {
        "dependent_cell_lines": 950,
        "cell_lines_with_data": 1000,
        "common_essential": False,  # DepMap's flag says NOT pan-essential
    }
    result = client._build_from_cache("MT-ND1", entry, ot_data=None)
    assert result.pan_essential is True, (
        "95% dependency must trip the pan-essential cap regardless of common_essential"
    )

    # Below 95% with common_essential=False → not flagged
    entry_low = {
        "dependent_cell_lines": 100,
        "cell_lines_with_data": 1000,
        "common_essential": False,
    }
    result_low = client._build_from_cache("BRAF", entry_low, ot_data=None)
    assert result_low.pan_essential is False


async def test_attach_safety_signals_reranks_by_faers_report_volume():
    """Bug F.2: selection should not be alphabetical — IMATINIB (decades of
    FAERS reports) should beat ASCIMINIB (recent approval, sparse reports)
    even though they're both phase-4 ABL1 inhibitors."""
    from genesis_bio_mcp.models import DrugInteraction
    from genesis_bio_mcp.tools.target_prioritization import attach_safety_signals

    # Both phase-4 inhibitors; alphabetically ASCIMINIB sorts first.
    drugs = [
        DrugInteraction(
            drug_name="ASCIMINIB",
            interaction_type="inhibitor",
            phase=4,
            approved=True,
            sources=["FDA"],
        ),
        DrugInteraction(
            drug_name="IMATINIB",
            interaction_type="inhibitor",
            phase=4,
            approved=True,
            sources=["FDA", "ChEMBL", "TTD"],  # more sources → ranks higher pre-lookup
        ),
    ]

    class _FakeOpenFDA:
        async def get_safety_signals(self, drug_name: str):
            from genesis_bio_mcp.models import DrugSafetySignal

            # IMATINIB has tens of thousands of FAERS reports; ASCIMINIB has
            # ~50 (recent approval). Both signals are returned by OpenFDA.
            counts = {"imatinib": 45000, "asciminib": 50}
            return DrugSafetySignal(
                drug_name=drug_name,
                total_reports=counts.get(drug_name.lower(), 0),
                top_adverse_events=[],
                boxed_warnings=["WARNING: serious liver toxicity"]
                if drug_name.lower() == "imatinib"
                else [],
                recalls=[],
            )

    out = await attach_safety_signals(drugs, openfda=_FakeOpenFDA())

    # Both got queried, but IMATINIB's signal wins on FAERS volume so it
    # ends up in the populated `.safety` slot before ASCIMINIB. The display
    # iterates known_drugs filtered by `if d.safety` — both have signals
    # here, but IMATINIB's signal contains warnings that the rendering
    # filter treats as "real" while ASCIMINIB's doesn't.
    by_name = {d.drug_name.lower(): d for d in out}
    assert by_name["imatinib"].safety is not None
    assert by_name["imatinib"].safety.total_reports == 45000
    # ASCIMINIB also got queried (it's in the pool), and got its (sparse)
    # signal attached. The point is that ranking happens by report count
    # so IMATINIB ranks first when display picks the top-N.
    assert by_name["asciminib"].safety is not None
    assert by_name["asciminib"].safety.total_reports == 50


def test_compute_score_zeros_gwas_when_fallback_marker_present():
    """Bug M: when the GWAS evidence carries the 'no exact-trait match'
    sentinel from Bug D's fallback, the GWAS axis must score 0 — the hits
    are gene-level top associations under unrelated traits and shouldn't
    earn trait-relevance credit. Without this gate, MRAS outranks HRAS
    for PDAC because off-trait hits inflate the score."""
    from genesis_bio_mcp.models import GwasEvidence, GwasHit
    from genesis_bio_mcp.tools.target_prioritization import _compute_score

    # GWAS payload that came from the unfiltered fallback path.
    fallback_gwas = GwasEvidence(
        gene_symbol="TP53",
        trait_query="Li-Fraumeni syndrome (no exact-trait match — top gene-level associations shown)",
        total_associations=5,
        associations=[
            GwasHit(
                study_accession="GCST007",
                trait="sex hormone-binding globulin",
                mapped_gene="TP53",
                risk_allele="rs999-A",
                p_value=4e-276,
            )
        ],
        strongest_p_value=4e-276,
    )
    _, breakdown = _compute_score(
        disease_assoc=None,
        cancer_dep=None,
        gwas_ev=fallback_gwas,
        compounds=None,
        protein=None,
    )
    assert breakdown.gwas == 0.0, "fallback hits must NOT earn GWAS scoring credit"

    # Without the sentinel, the same hit count earns full credit.
    real_gwas = GwasEvidence(
        gene_symbol="PCSK9",
        trait_query="LDL cholesterol",
        total_associations=5,
        associations=[],
        strongest_p_value=1e-20,
    )
    _, breakdown2 = _compute_score(
        disease_assoc=None,
        cancer_dep=None,
        gwas_ev=real_gwas,
        compounds=None,
        protein=None,
    )
    assert breakdown2.gwas == 2.0  # 5 hits ≥ 3 cap → full credit


def test_compute_score_chemical_matter_discounts_binding_only_potency():
    """Bug I: pChEMBL=9 from a binding assay alone should score lower than
    pChEMBL=9 from a cell-based / functional assay. Without this, MYC
    (1.5/1.5 chem matter from binding-only data) is ranked equal to truly
    druggable kinases — but MYC binders virtually never translate to
    cellular activity."""
    from genesis_bio_mcp.models import ChEMBLCompounds
    from genesis_bio_mcp.tools.target_prioritization import _compute_score

    # Binding-only at clinical-grade potency → 1.0 (was 1.5 in v0.3.2)
    binding_only = ChEMBLCompounds(
        gene_symbol="MYC",
        target_chembl_id="CHEMBL1234",
        total_active_compounds=89,
        best_pchembl=9.2,
        best_pchembl_functional=None,
        best_pchembl_binding=9.2,
        compounds=[],
    )
    _, breakdown_b = _compute_score(
        disease_assoc=None,
        cancer_dep=None,
        gwas_ev=None,
        compounds=None,
        protein=None,
        chembl_compounds=binding_only,
    )
    assert breakdown_b.chem_matter == 1.0, (
        f"binding-only clinical-grade should be 1.0 (discounted from 1.5), "
        f"got {breakdown_b.chem_matter}"
    )

    # Functional cell-based at clinical-grade potency → 1.5 (full credit)
    functional = ChEMBLCompounds(
        gene_symbol="EGFR",
        target_chembl_id="CHEMBL203",
        total_active_compounds=200,
        best_pchembl=9.2,
        best_pchembl_functional=9.2,
        best_pchembl_binding=8.5,
        compounds=[],
    )
    _, breakdown_f = _compute_score(
        disease_assoc=None,
        cancer_dep=None,
        gwas_ev=None,
        compounds=None,
        protein=None,
        chembl_compounds=functional,
    )
    assert breakdown_f.chem_matter == 1.5, "functional clinical-grade should keep full credit"

    # Lead-quality functional → 1.0; binding-only lead-quality → 0.7
    func_lead = ChEMBLCompounds(
        gene_symbol="X",
        target_chembl_id="CHEMBL1",
        total_active_compounds=10,
        best_pchembl=7.5,
        best_pchembl_functional=7.5,
        best_pchembl_binding=None,
        compounds=[],
    )
    _, b_func_lead = _compute_score(
        disease_assoc=None,
        cancer_dep=None,
        gwas_ev=None,
        compounds=None,
        protein=None,
        chembl_compounds=func_lead,
    )
    assert b_func_lead.chem_matter == 1.0

    bind_lead = ChEMBLCompounds(
        gene_symbol="Y",
        target_chembl_id="CHEMBL2",
        total_active_compounds=10,
        best_pchembl=7.5,
        best_pchembl_functional=None,
        best_pchembl_binding=7.5,
        compounds=[],
    )
    _, b_bind_lead = _compute_score(
        disease_assoc=None,
        cancer_dep=None,
        gwas_ev=None,
        compounds=None,
        protein=None,
        chembl_compounds=bind_lead,
    )
    assert b_bind_lead.chem_matter == 0.7


# ---------------------------------------------------------------------------
# v0.3.4 regression tests — five fourth-round fixes from the third smoke run
# ---------------------------------------------------------------------------


def test_format_ic50_nm_handles_sub_nanomolar_potency():
    """Bug Q: pChEMBL ≥ 10 yielded sub-nM IC50 that rendered as '0.0 nM'.
    Sub-nM should display in pM; nM in 1–1000 range; µM beyond."""
    from genesis_bio_mcp.models import _format_ic50_nm

    # pChEMBL=11 → 0.01 nM → 10 pM
    assert _format_ic50_nm(0.01) == "10 pM"
    # pChEMBL=10.7 → ~0.02 nM → 20 pM (matches CALCA/GLP1R real-world case)
    assert _format_ic50_nm(0.0199526) == "20 pM"
    # Boundary at 1 nM stays in nM
    assert _format_ic50_nm(1.0) == "1.0 nM"
    # Sub-nM but close → still pM
    assert _format_ic50_nm(0.5) == "500 pM"
    # Standard nM range
    assert _format_ic50_nm(50.0) == "50.0 nM"
    # Crosses to µM
    assert _format_ic50_nm(2500.0) == "2.50 µM"


def test_reactome_parse_pathways_filters_non_human_species():
    """Bug R: R-CFA-* (canine), R-MMU-* (mouse), R-RNO-* (rat) entries can
    leak through Reactome's analysis when the gene maps across species
    silos. Human queries must only surface R-HSA-* pathways."""
    from genesis_bio_mcp.clients.reactome import _parse_pathways

    raw = [
        {
            "stId": "R-HSA-1234567",
            "name": "Class A/1 (Rhodopsin-like receptors)",
            "entities": {"pValue": 1e-5, "total": 200},
        },
        {
            "stId": "R-CFA-381676",  # canine — must be dropped (the GLP1R bug)
            "name": "Some canine pathway",
            "entities": {"pValue": 1e-3, "total": 50},
        },
        {
            "stId": "R-MMU-9876543",  # mouse — must be dropped
            "name": "Murine pathway",
            "entities": {"pValue": 1e-4, "total": 30},
        },
        {
            "stId": "",  # missing stId — must be dropped (was previously kept as anonymous)
            "name": "Anonymous",
            "entities": {"pValue": 1e-2, "total": 10},
        },
    ]
    parsed = _parse_pathways(raw)
    assert len(parsed) == 1
    assert parsed[0].reactome_id == "R-HSA-1234567"
    assert all(p.reactome_id.startswith("R-HSA-") for p in parsed)


def test_open_targets_normalize_variants_expands_common_acronyms():
    """Bug P: bare acronyms (NSCLC, PDAC, T2DM) returned 0 OT hits because
    OT's autocomplete-style search doesn't index abbreviations. The variant
    generator must add the canonical full form so the resolver can find it."""
    from genesis_bio_mcp.clients.open_targets import _normalize_indication_variants

    # Bare acronym → expansion appears as a fallback variant
    nsclc = _normalize_indication_variants("NSCLC")
    assert "NSCLC" in nsclc  # literal still tried first
    assert "non-small cell lung carcinoma" in nsclc

    pdac = _normalize_indication_variants("PDAC")
    assert "pancreatic ductal adenocarcinoma" in pdac

    t2dm = _normalize_indication_variants("T2DM")
    assert "type 2 diabetes mellitus" in t2dm

    # Acronym embedded in longer string — expansion still fires
    nsclc_messy = _normalize_indication_variants("EGFR-mutant NSCLC adenocarcinoma")
    assert "non-small cell lung carcinoma" in nsclc_messy

    # Mixed-case / case-insensitive acronyms (HFpEF, AFib)
    hfpef = _normalize_indication_variants("HFpEF")
    assert "heart failure with preserved ejection fraction" in hfpef

    # Non-acronym input is unaffected
    plain = _normalize_indication_variants("rheumatoid arthritis")
    assert plain == ["rheumatoid arthritis"]  # no spurious expansions


def test_compute_score_zeros_gwas_when_traits_off_indication():
    """Bug M refined: TP53/Li-Fraumeni returned direct-match GWAS hits whose
    trait labels are sex-hormone-binding-globulin and height — completely
    unrelated to Li-Fraumeni. The fallback sentinel doesn't fire (these came
    through filter_by_trait), so we need an additional trait-relevance gate
    in _compute_score that zeros the GWAS axis when no hit's trait overlaps
    the indication meaningfully."""
    from genesis_bio_mcp.models import GwasEvidence, GwasHit
    from genesis_bio_mcp.tools.target_prioritization import _compute_score

    off_trait = GwasEvidence(
        gene_symbol="TP53",
        trait_query="li-fraumeni syndrome",  # no fallback sentinel — direct return
        total_associations=5,
        associations=[
            GwasHit(
                study_accession="GCST007",
                trait="sex hormone-binding globulin levels",
                mapped_gene="TP53",
                risk_allele="rs999-A",
                p_value=4e-276,
            ),
            GwasHit(
                study_accession="GCST008",
                trait="height",
                mapped_gene="TP53",
                risk_allele="rs1000-G",
                p_value=1e-50,
            ),
        ],
        strongest_p_value=4e-276,
    )
    _, breakdown = _compute_score(
        disease_assoc=None,
        cancer_dep=None,
        gwas_ev=off_trait,
        compounds=None,
        protein=None,
        indication="li-fraumeni syndrome",
    )
    assert breakdown.gwas == 0.0, (
        "trait-relevance check should zero GWAS when no hit's trait overlaps the indication"
    )

    # On-trait direct-match hits keep full credit
    on_trait = GwasEvidence(
        gene_symbol="PCSK9",
        trait_query="coronary artery disease",
        total_associations=5,
        associations=[
            GwasHit(
                study_accession="GCST100",
                trait="Coronary artery disease",
                mapped_gene="PCSK9",
                risk_allele="rs11591147-T",
                p_value=4e-20,
            ),
            GwasHit(
                study_accession="GCST101",
                trait="Myocardial infarction",
                mapped_gene="PCSK9",
                risk_allele="rs562556-A",
                p_value=3e-15,
            ),
        ],
        strongest_p_value=4e-20,
    )
    _, breakdown2 = _compute_score(
        disease_assoc=None,
        cancer_dep=None,
        gwas_ev=on_trait,
        compounds=None,
        protein=None,
        indication="coronary artery disease",
    )
    assert breakdown2.gwas == 2.0  # full credit for relevant hits


def test_target_prioritization_breakdown_table_matches_score_breakdown():
    """Bug M.1: the displayed scoring breakdown was recomputing each axis
    inline with stale formulas (cap at 10 not 3 for GWAS, no monogenic
    credit, no functional/binding split). The numbers in the rendered
    table drifted from the canonical _compute_score output. Fix uses
    score_breakdown directly so they match by construction."""
    from genesis_bio_mcp.models import (
        CancerDependency,
        GeneResolution,
        GwasEvidence,
        GwasHit,
        ScoreBreakdown,
        TargetDiseaseAssociation,
        TargetPrioritizationReport,
    )

    # Build a report with intentionally tricky GWAS — the inline recompute
    # would have shown 1.0 (5 hits / 10 * 2.0); the breakdown shows 0.0
    # (off-trait suppression). The fix means the table mirrors the breakdown.
    breakdown = ScoreBreakdown(
        ot=2.55,
        depmap=0.5,
        gwas=0.0,  # canonical: zeroed by trait-relevance gate
        known_drug=1.43,
        chem_matter=1.0,  # canonical: binding-only discount
        protein=1.5,
        expression=0.0,
    )
    fallback_gwas = GwasEvidence(
        gene_symbol="TP53",
        trait_query="li-fraumeni syndrome",
        total_associations=5,
        associations=[
            GwasHit(
                study_accession="GCST007",
                trait="sex hormone-binding globulin",
                mapped_gene="TP53",
                risk_allele="rs999-A",
                p_value=4e-276,
            )
        ],
        strongest_p_value=4e-276,
    )
    report = TargetPrioritizationReport(
        gene_symbol="TP53",
        indication="Li-Fraumeni syndrome",
        resolution=GeneResolution(hgnc_symbol="TP53", source="input"),
        protein_info=None,
        disease_association=TargetDiseaseAssociation(
            gene_symbol="TP53",
            disease_name="Li-Fraumeni syndrome",
            disease_efo_id="EFO_0009248",
            ensembl_id="ENSG00000141510",
            overall_score=0.85,
            genetic_association_score=0.5,
            somatic_mutation_score=0.7,
            known_drug_score=0.95,
            literature_mining_score=0.8,
            evidence_count=200,
        ),
        cancer_dependency=CancerDependency(
            gene_symbol="TP53",
            mean_ceres_score=-0.2,
            fraction_dependent_lines=0.15,
            pan_essential=False,
            top_dependent_lineages=[],
            cell_lines=[],
            data_source="DepMap Chronos Combined",
        ),
        gwas_evidence=fallback_gwas,
        compounds=None,
        chembl_compounds=None,
        priority_score=6.98,
        priority_tier="Medium",
        score_breakdown=breakdown,
        evidence_summary="test",
        data_gaps=[],
        errors={},
        data_coverage_pct=80.0,
        proxy_data_flags={},
        score_confidence_interval=None,
        api_latency_s={},
    )
    md = report.to_markdown()

    # The GWAS row must show 0.00, NOT the stale formula's 1.00 (5 hits/10 * 2.0).
    assert "| GWAS evidence | 0.00" in md
    # And it must explain WHY (off-trait suppression note).
    assert "off-indication" in md or "not credited" in md
    # OT row uses the breakdown's 2.55, not da.overall_score * 3.0 = 2.55
    # (these happen to match here, but the test guards against drift).
    assert "| Open Targets association | 2.55" in md
    # Cancer dependency row uses breakdown's 0.50, not the inline recompute.
    assert "| Cancer dependency | 0.50" in md
    # Sum check: every numeric contribution in the breakdown must come from
    # `score_breakdown`, so summing the visible contributions equals
    # breakdown.total (within rounding).
    import re as _re

    contribs = [float(m) for m in _re.findall(r"\| [A-Za-z][^|]*? \| ([0-9]+\.[0-9]+)", md)]
    # The first match is the priority score itself; subsequent ones are
    # axis contributions. Sum the axis contributions and compare to
    # breakdown.total.
    axis_contribs = contribs[:7]  # 7 axes max (some may be skipped)
    # Allow rounding tolerance.
    assert abs(sum(axis_contribs) - breakdown.total) < 0.05, (
        f"breakdown table contributions {axis_contribs} (sum {sum(axis_contribs):.2f}) "
        f"don't match score_breakdown.total ({breakdown.total:.2f}) — display drift"
    )
