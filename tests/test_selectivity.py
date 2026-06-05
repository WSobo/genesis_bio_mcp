"""Tests for assess_selectivity / the ChEMBL off-target lookup.

ChEMBL HTTP is mocked with respx so these are deterministic and CI-safe (ChEMBL itself is
flaky under load). The corpus predicted layer is exercised separately in test_corpus.py.
"""

from __future__ import annotations

import httpx
import respx

from genesis_bio_mcp.clients.chembl import ChEMBLClient
from genesis_bio_mcp.tools.selectivity import assess_selectivity

# A real, promiscuous kinase inhibitor — parses fine in RDKit (we only mock the HTTP).
DASATINIB = "Cc1nc(Nc2ncc(C(=O)Nc3c(C)cccc3Cl)s2)cc(N2CCN(CCO)CC2)n1"
_MOL_RE = r"chembl/api/data/molecule/[\w-]+\.json"
_ACT = "https://www.ebi.ac.uk/chembl/api/data/activity"


def _activity_payload(targets: list[tuple[str, str, float]]) -> dict:
    """Build a ChEMBL /activity response from (target_chembl_id, pref_name, pchembl) rows."""
    return {
        "activities": [
            {
                "target_chembl_id": tid,
                "target_pref_name": name,
                "target_organism": "Homo sapiens",
                "pchembl_value": str(p),
            }
            for tid, name, p in targets
        ],
        "page_meta": {"total_count": len(targets)},
    }


@respx.mock
async def test_promiscuous_compound_profile():
    respx.get(url__regex=_MOL_RE).mock(
        return_value=httpx.Response(200, json={"molecule_chembl_id": "CHEMBL1421"})
    )
    targets = [
        ("CHEMBL1862", "Tyrosine-protein kinase ABL1", 10.0),
        ("CHEMBL267", "Proto-oncogene tyrosine-protein kinase Src", 9.1),
        ("CHEMBL2068", "Macrophage colony-stimulating factor 1 receptor", 8.7),
        ("CHEMBL258", "Tyrosine-protein kinase Lyn", 8.6),
        ("CHEMBL1841", "Ephrin type-A receptor 3", 8.5),
        ("CHEMBL1936", "Stem cell growth factor receptor KIT", 8.4),
    ]
    respx.get(url__startswith=_ACT).mock(
        return_value=httpx.Response(200, json=_activity_payload(targets))
    )

    async with httpx.AsyncClient() as c:
        prof = await assess_selectivity(DASATINIB, chembl=ChEMBLClient(c), corpus_pool=None)

    assert prof is not None
    assert prof.molecule_chembl_id == "CHEMBL1421"
    assert len(prof.measured_off_targets) == 6
    assert prof.measured_off_targets[0].best_pchembl == 10.0  # sorted by potency
    assert prof.selectivity_label == "Promiscuous"  # 5 potent (>=7) off-targets
    assert prof.corpus_used is False
    md = prof.to_markdown()
    assert "ABL1" in md and "Src" in md
    assert "absence ≠ selectivity" in md  # honesty caveat


@respx.mock
async def test_selective_compound_profile():
    respx.get(url__regex=_MOL_RE).mock(
        return_value=httpx.Response(200, json={"molecule_chembl_id": "CHEMBLX"})
    )
    respx.get(url__startswith=_ACT).mock(
        return_value=httpx.Response(
            200, json=_activity_payload([("CHEMBL1", "Primary kinase", 9.0)])
        )
    )
    async with httpx.AsyncClient() as c:
        prof = await assess_selectivity("c1ccccc1", chembl=ChEMBLClient(c), corpus_pool=None)
    assert len(prof.measured_off_targets) == 1
    assert prof.selectivity_label == "Selective"
    assert "Primary kinase" in (prof.primary_target or "")


@respx.mock
async def test_compound_not_in_chembl():
    respx.get(url__regex=_MOL_RE).mock(return_value=httpx.Response(404))
    async with httpx.AsyncClient() as c:
        prof = await assess_selectivity("c1ccccc1", chembl=ChEMBLClient(c), corpus_pool=None)
    assert prof.molecule_chembl_id is None
    assert prof.measured_off_targets == []
    assert prof.selectivity_label == "Unknown (no measured data)"
    assert "not found in ChEMBL" in prof.to_markdown()


async def test_invalid_smiles_returns_none():
    async with httpx.AsyncClient() as c:
        assert (
            await assess_selectivity("not_a_smiles", chembl=ChEMBLClient(c), corpus_pool=None)
            is None
        )


@respx.mock
async def test_activity_500_retries_then_recovers(monkeypatch):
    """A transient ChEMBL 5xx is retried, not fatal."""
    import genesis_bio_mcp.clients.chembl as ch

    async def _no_sleep(*_a, **_k):
        return None

    monkeypatch.setattr(ch.asyncio, "sleep", _no_sleep)
    respx.get(url__regex=_MOL_RE).mock(
        return_value=httpx.Response(200, json={"molecule_chembl_id": "CHEMBL1421"})
    )
    respx.get(url__startswith=_ACT).mock(
        side_effect=[
            httpx.Response(500),
            httpx.Response(200, json=_activity_payload([("CHEMBL1", "ABL1 kinase", 9.0)])),
        ]
    )
    async with httpx.AsyncClient() as c:
        prof = await assess_selectivity(DASATINIB, chembl=ChEMBLClient(c), corpus_pool=None)
    assert prof.molecule_chembl_id == "CHEMBL1421"
    assert len(prof.measured_off_targets) == 1
