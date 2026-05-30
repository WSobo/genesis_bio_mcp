"""End-to-end tests for the target prioritization orchestration tool.

These tests mock at the client level (no real HTTP) to validate the full
data flow from tool invocation to structured TargetPrioritizationReport.
"""

from unittest.mock import AsyncMock, MagicMock

import pytest

from genesis_bio_mcp.tools.gene_resolver import resolve_gene
from genesis_bio_mcp.tools.target_prioritization import prioritize_target
from tests.conftest import (
    build_mock_association,
    build_mock_compounds,
    build_mock_dependency,
    build_mock_gwas,
    build_mock_protein_info,
)


def _mock_clients(
    protein=None,
    association=None,
    dependency=None,
    gwas=None,
    compounds=None,
    resolution=None,
):
    """Build a set of AsyncMock clients with pre-configured return values."""
    uniprot = AsyncMock()
    uniprot.get_protein.return_value = protein or build_mock_protein_info()
    # _search is used by gene_resolver internally
    uniprot._search.return_value = {
        "primaryAccession": "P15056",
        "genes": [
            {
                "geneName": {"value": "BRAF"},
                "synonyms": [{"value": "BRAF1"}],
            }
        ],
        "entryType": "UniProtKB reviewed (Swiss-Prot)",
    }
    uniprot.search_by_synonym.return_value = None
    # Use MagicMock for HTTP responses so raise_for_status() stays synchronous
    _ncbi_response = MagicMock()
    _ncbi_response.json.return_value = {"esearchresult": {"idlist": ["673"]}}
    _ncbi_response.raise_for_status = MagicMock()
    _ncbi_client = AsyncMock()
    _ncbi_client.get.return_value = _ncbi_response
    uniprot._client = _ncbi_client

    open_targets = AsyncMock()
    open_targets.get_association.return_value = association or build_mock_association()

    depmap = AsyncMock()
    depmap.get_essentiality.return_value = dependency or build_mock_dependency()

    gwas_client = AsyncMock()
    gwas_client.get_evidence.return_value = gwas or build_mock_gwas()

    pubchem = AsyncMock()
    pubchem.get_compounds.return_value = compounds or build_mock_compounds()

    chembl = AsyncMock()
    chembl.get_compounds.return_value = None  # ChEMBL absent by default; tests opt-in

    return uniprot, open_targets, depmap, gwas_client, pubchem, chembl


# ---------------------------------------------------------------------------
# Happy path: BRAF / melanoma full report
# ---------------------------------------------------------------------------


async def test_braf_melanoma_full_report():
    """BRAF assessed for melanoma should return a strong evidence profile with all fields."""
    uniprot, open_targets, depmap, gwas_client, pubchem, chembl = _mock_clients(
        protein=build_mock_protein_info("BRAF"),
        association=build_mock_association("BRAF", "melanoma", score=0.89),
        dependency=build_mock_dependency("BRAF", mean_score=-0.72, fraction_dependent=0.61),
        gwas=build_mock_gwas("BRAF", "melanoma", n_hits=5),
        compounds=build_mock_compounds("BRAF", n_active=312),
    )

    report = await prioritize_target(
        gene_symbol="BRAF",
        indication="melanoma",
        uniprot=uniprot,
        open_targets=open_targets,
        depmap=depmap,
        gwas=gwas_client,
        pubchem=pubchem,
        chembl=chembl,
    )

    # Core identity
    assert report.gene_symbol == "BRAF"
    assert report.indication == "melanoma"

    # Evidence profile — each axis rated independently; no composite score
    assert not hasattr(report, "priority_score")
    axes = {a.name: a.rating for a in report.evidence_profile}
    assert axes["Disease association"] == "strong"  # OT overall 0.89
    assert axes["Clinical validation"] == "strong"  # known-drug 0.88
    assert axes["Genetic evidence"] == "strong"  # OT genetic 0.72 + 5 melanoma GWAS hits
    assert axes["Cancer dependency"] == "strong"  # 61% dependent, real DepMap Chronos

    # Evidence fields all populated
    assert report.protein_info is not None
    assert report.protein_info.reviewed is True
    assert report.disease_association is not None
    assert report.disease_association.overall_score == pytest.approx(0.89)
    assert report.cancer_dependency is not None
    assert report.cancer_dependency.fraction_dependent_lines == pytest.approx(0.61)
    assert report.gwas_evidence is not None
    assert report.gwas_evidence.total_associations == 5
    assert report.compounds is not None
    assert report.compounds.total_active_compounds == 312

    # Core APIs all succeeded — only ChEMBL (optional, returns None by default in mock) may be a gap
    assert set(report.data_gaps) <= {"chembl"}
    assert len(report.errors) == 0

    # Evidence summary is non-empty and mentions BRAF
    assert len(report.evidence_summary) > 50
    assert "BRAF" in report.evidence_summary
    assert "melanoma" in report.evidence_summary.lower()

    # Resolution populated
    assert report.resolution.hgnc_symbol == "BRAF"


# ---------------------------------------------------------------------------
# Graceful partial failure
# ---------------------------------------------------------------------------


async def test_prioritization_handles_api_failures_gracefully():
    """When DepMap and GWAS fail, report still returns with correct data_gaps."""
    uniprot, open_targets, depmap, gwas_client, pubchem, chembl = _mock_clients()

    depmap.get_essentiality.side_effect = Exception("Connection timeout")
    gwas_client.get_evidence.side_effect = Exception("503 Service Unavailable")

    report = await prioritize_target(
        gene_symbol="BRAF",
        indication="melanoma",
        uniprot=uniprot,
        open_targets=open_targets,
        depmap=depmap,
        gwas=gwas_client,
        pubchem=pubchem,
        chembl=chembl,
    )

    # Failed modules result in None fields
    assert report.cancer_dependency is None
    assert report.gwas_evidence is None

    # data_gaps tracks which modules had no data
    assert "depmap" in report.data_gaps
    assert "gwas" in report.data_gaps

    # Errors recorded
    assert "depmap" in report.errors
    assert "gwas" in report.errors

    # Remaining evidence still present
    assert report.disease_association is not None
    assert report.compounds is not None

    # Failed sources are rated n/a, not silently scored zero
    axes = {a.name: a.rating for a in report.evidence_profile}
    assert axes["Cancer dependency"] == "n/a"  # DepMap failed
    assert axes["Genetic evidence"] == "strong"  # GWAS failed but OT genetic 0.72 still rates it
    assert axes["Disease association"] == "strong"  # OT survived

    # Does not raise — critical for agent stability
    assert report is not None


async def test_prioritization_all_apis_fail():
    """When all APIs fail, report returns with all data_gaps and Low priority."""
    uniprot, open_targets, depmap, gwas_client, pubchem, chembl = _mock_clients()

    uniprot.get_protein.side_effect = Exception("UniProt unreachable")
    open_targets.get_association.side_effect = Exception("Open Targets unreachable")
    depmap.get_essentiality.side_effect = Exception("DepMap unreachable")
    gwas_client.get_evidence.side_effect = Exception("GWAS unreachable")
    pubchem.get_compounds.side_effect = Exception("PubChem unreachable")
    chembl.get_compounds.side_effect = Exception("ChEMBL unreachable")

    report = await prioritize_target(
        gene_symbol="TP53",
        indication="lung cancer",
        uniprot=uniprot,
        open_targets=open_targets,
        depmap=depmap,
        gwas=gwas_client,
        pubchem=pubchem,
        chembl=chembl,
    )

    assert report.protein_info is None
    assert report.disease_association is None
    assert report.cancer_dependency is None
    assert report.gwas_evidence is None
    assert report.compounds is None
    assert len(report.data_gaps) == 6  # includes chembl
    # Every axis is n/a — no data means no signal, never a fabricated zero score
    assert report.evidence_profile
    assert all(a.rating == "n/a" for a in report.evidence_profile)
    assert report is not None  # Never raises


# ---------------------------------------------------------------------------
# Score computation edge cases
# ---------------------------------------------------------------------------


async def test_pan_essential_gene_rated_weak_not_strong():
    """Pan-essential genes carry a narrow therapeutic window — the cancer-
    dependency axis must be 'weak' (a liability), never 'strong', no matter
    how high the dependent fraction is."""
    pan_essential_dep = build_mock_dependency("MYC", mean_score=-1.8, fraction_dependent=0.95)
    pan_essential_dep = pan_essential_dep.model_copy(update={"pan_essential": True})

    uniprot, open_targets, depmap, gwas_client, pubchem, chembl = _mock_clients(
        dependency=pan_essential_dep
    )

    report = await prioritize_target(
        gene_symbol="MYC",
        indication="lymphoma",
        uniprot=uniprot,
        open_targets=open_targets,
        depmap=depmap,
        gwas=gwas_client,
        pubchem=pubchem,
        chembl=chembl,
    )

    assert report is not None
    dep_axis = next(a for a in report.evidence_profile if a.name == "Cancer dependency")
    assert dep_axis.rating == "weak"  # pan-essential is a liability, not strength
    assert "pan-essential" in dep_axis.detail.lower()


# ---------------------------------------------------------------------------
# Gene resolution integration
# ---------------------------------------------------------------------------


async def test_gene_resolver_uses_canonical_symbol():
    """resolve_gene should return the canonical HGNC symbol, not the input alias."""
    uniprot_mock = AsyncMock()
    uniprot_mock._search.return_value = {
        "primaryAccession": "P15056",
        "genes": [
            {
                "geneName": {"value": "BRAF"},
                "synonyms": [{"value": "BRAF1"}],
            }
        ],
    }
    uniprot_mock.search_by_synonym.return_value = None
    # Use MagicMock for the HTTP response so raise_for_status() stays synchronous
    _ncbi_response = MagicMock()
    _ncbi_response.json.return_value = {"esearchresult": {"idlist": ["673"]}}
    _ncbi_response.raise_for_status = MagicMock()
    _ncbi_client = AsyncMock()
    _ncbi_client.get.return_value = _ncbi_response
    uniprot_mock._client = _ncbi_client

    result = await resolve_gene("braf1", uniprot_client=uniprot_mock)

    assert result.hgnc_symbol == "BRAF"
    assert result.uniprot_accession == "P15056"
