"""Regression tests: response parsing tolerates sparse/malformed list entries.

Several parsers used hard dict subscripts (e.g. ``d["id"]``) on API-response
fields that sit *outside* the HTTP try/except, so a valid-but-sparse response
(a list entry missing an expected key) raised an unhandled ``KeyError`` that
escaped to the MCP tool caller. The fixes switch to guarded ``.get()``; these
tests assert a malformed entry is skipped rather than fatal.
"""

from __future__ import annotations

import httpx
import respx

from genesis_bio_mcp.clients.interpro import InterProClient
from genesis_bio_mcp.clients.open_targets import _parse_row


def test_parse_row_tolerates_datatypescore_missing_id() -> None:
    """An Open Targets datatypeScores entry without 'id' is skipped, not fatal."""
    row = {
        "score": 0.7,
        "datatypeScores": [
            {"id": "genetic_association", "score": 0.5},
            {"score": 0.3},  # malformed: no "id" — used to raise KeyError
        ],
    }
    assoc = _parse_row(row, "EGFR", "lung cancer", "EFO_0000001", "ENSG00000146648")
    assert assoc.overall_score == 0.7
    assert assoc.genetic_association_score == 0.5
    assert assoc.evidence_count == 1  # only the well-formed entry survives


@respx.mock
async def test_interpro_tolerates_go_term_missing_name() -> None:
    """An InterPro GO term missing 'name' renders without crashing."""
    payload = {
        "count": 1,
        "results": [
            {
                "metadata": {
                    "accession": "IPR000001",
                    "name": "Test domain",
                    "type": "domain",
                    "go_terms": [{"identifier": "GO:0005524"}],  # missing "name"
                },
                "proteins": [
                    {"entry_protein_locations": [{"fragments": [{"start": 10, "end": 50}]}]}
                ],
            }
        ],
    }
    respx.get(url__regex=r".*interpro.*").mock(return_value=httpx.Response(200, json=payload))
    async with httpx.AsyncClient() as http:
        result = await InterProClient(http).get_domains("BRAF", "P15056")
    assert result is not None
    assert result.total_entries == 1
    assert len(result.domains) == 1
    assert result.domains[0].go_terms == ["GO:0005524"]
