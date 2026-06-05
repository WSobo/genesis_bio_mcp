"""Regression tests for data-quality fixes found during live testing.

1. UniProt ``get_protein`` returned a *moonlighting* function (e.g. EGFR → an
   HCV-receptor note) because the parser kept the LAST FUNCTION comment; it must
   keep the first (primary) one.
2. GTEx silently dropped tissues whose median TPM is exactly ``0.0`` because of a
   ``median or medianTpm`` falsy fallthrough; 0.0 must be retained.
"""

from __future__ import annotations

from unittest.mock import MagicMock

import httpx
import respx

from genesis_bio_mcp.clients.gtex import GTExClient
from genesis_bio_mcp.clients.uniprot import _parse_entry


def test_uniprot_function_summary_takes_first_function_comment() -> None:
    entry = {
        "primaryAccession": "P00533",
        "entryType": "UniProtKB reviewed (Swiss-Prot)",
        "comments": [
            {
                "commentType": "FUNCTION",
                "texts": [{"value": "Receptor tyrosine kinase binding EGF."}],
            },
            {
                "commentType": "FUNCTION",
                "texts": [{"value": "(Microbial infection) Acts as a receptor for HCV."}],
            },
        ],
    }
    info = _parse_entry(entry, "EGFR")
    assert info.function_summary == "Receptor tyrosine kinase binding EGF."


@respx.mock
async def test_gtex_keeps_zero_median_tissues() -> None:
    payload = {
        "data": [
            {"tissueSiteDetailId": "Pancreas", "median": 100.0, "sampleCount": 5},
            {"tissueSiteDetailId": "Whole_Blood", "median": 0.0, "sampleCount": 5},
        ]
    }
    respx.get(url__regex=r".*medianGeneExpression.*").mock(
        return_value=httpx.Response(200, json=payload)
    )
    async with httpx.AsyncClient() as http:
        client = GTExClient(http, ensembl=MagicMock())
        rows = await client._fetch_expression("ENSG00000000000.1")
    by_tissue = {r.tissue: r.median_tpm for r in rows}
    assert by_tissue["Whole - Blood"] == 0.0  # 0.0 retained, not dropped
    assert by_tissue["Pancreas"] == 100.0
