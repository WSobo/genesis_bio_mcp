"""Unit tests for the M5 provenance envelope (agent contract).

Every ``response_format="json"`` result is wrapped in a consistent
``{"provenance": {...}, "data": {...}}`` envelope (or ``{"provenance", "error"}``
on failure); Markdown output is unchanged.
"""

from __future__ import annotations

import json
from datetime import datetime
from types import SimpleNamespace

import pytest

from genesis_bio_mcp.models import GeneResolution, Provenance
from genesis_bio_mcp.server import (
    ComputeMolecularPropertiesInput,
    _build_provenance,
    _fmt,
    _now_iso,
)
from genesis_bio_mcp.server import (
    compute_molecular_properties as compute_molecular_properties_tool,
)
from genesis_bio_mcp.tools.cheminformatics import compute_molecular_properties


def test_json_success_envelope_wraps_data_and_provenance():
    props = compute_molecular_properties("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
    out = _fmt(props, "json", "")
    env = json.loads(out)
    assert set(env.keys()) == {"provenance", "data"}
    assert env["data"]["molecular_formula"] == "C9H8O4"
    prov = env["provenance"]
    assert prov["source"] == "RDKit (local)"
    assert prov["query"] == "CC(=O)Oc1ccccc1C(=O)O"  # pulled from input_smiles
    # retrieved_at is an ISO-8601 'Z' UTC timestamp
    assert prov["retrieved_at"].endswith("Z")
    datetime.strptime(prov["retrieved_at"], "%Y-%m-%dT%H:%M:%SZ")


def test_json_error_envelope_carries_provenance_and_typed_status():
    out = _fmt(None, "json", "No data found for FOO.")
    env = json.loads(out)
    assert set(env.keys()) == {"provenance", "error"}
    # M6: error is a typed object, not a bare string.
    assert env["error"] == {"status": "NotFound", "message": "No data found for FOO."}
    # Source is unknown for a missing result → generic server label, never crashes.
    assert env["provenance"]["source"] == "genesis-bio-mcp"
    assert env["provenance"]["query"] is None


def test_json_error_status_defaults_to_not_found():
    env = json.loads(_fmt(None, "json", "nope"))
    assert env["error"]["status"] == "NotFound"


def test_json_error_status_can_be_overridden():
    for status in ("InvalidInput", "RateLimited", "UpstreamUnavailable"):
        env = json.loads(_fmt(None, "json", "boom", error_status=status))
        assert env["error"]["status"] == status
        assert env["error"]["message"] == "boom"


def test_markdown_error_unaffected_by_status():
    # The typed status is JSON-only; Markdown error text stays clean.
    assert _fmt(None, "markdown", "boom", error_status="InvalidInput") == "**Error:** boom"


def test_markdown_output_is_unchanged_by_provenance():
    props = compute_molecular_properties("CCO")
    assert _fmt(props, "markdown", "") == props.to_markdown()
    assert _fmt(None, "markdown", "boom") == "**Error:** boom"


def test_provenance_source_and_query_from_gene_model():
    prov = _build_provenance(GeneResolution(hgnc_symbol="ERBB2", source="uniprot"))
    assert isinstance(prov, Provenance)
    assert prov.source == "UniProt + NCBI"
    assert prov.query == "ERBB2"  # from hgnc_symbol


def test_provenance_unmapped_model_falls_back_gracefully():
    # An object whose type is not in the source map → generic label, no crash.
    prov = _build_provenance(SimpleNamespace(mean_confidence=0.87, release="2024Q2"))
    assert prov.source == "genesis-bio-mcp"
    assert prov.confidence == 0.87  # derived from mean_confidence
    assert prov.source_version == "2024Q2"  # derived from release


def test_now_iso_is_utc_z_format():
    ts = _now_iso()
    assert ts.endswith("Z")
    datetime.strptime(ts, "%Y-%m-%dT%H:%M:%SZ")


@pytest.mark.asyncio
async def test_tool_emits_invalid_input_status_for_bad_smiles():
    # End-to-end: a real tool classifies a malformed SMILES as InvalidInput, not NotFound.
    out = await compute_molecular_properties_tool(
        ComputeMolecularPropertiesInput(smiles="not_a_smiles!!!", response_format="json")
    )
    env = json.loads(out)
    assert env["error"]["status"] == "InvalidInput"
    assert "could not be parsed" in env["error"]["message"]
