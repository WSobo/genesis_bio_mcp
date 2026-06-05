"""Data-quality round 2 regression tests (found during live testing).

1. MaveDB ``get_dms_scores`` returned off-gene datasets because the text search
   matches the gene name in unrelated studies' abstracts; must filter to the
   queried target gene.
2. UniProt variant parsing emitted a bare "30" for features with no substitution
   and labelled every non-disease polymorphism "unknown disease"; must skip empty
   features, extract the disease token, and not invent diseases.
3. ``search_similar_compounds`` claimed "Tanimoto >= threshold" while the displayed
   column is a different (RDKit ECFP4) metric; the header must say so honestly.
"""

from __future__ import annotations

import httpx
import respx

from genesis_bio_mcp.clients.mavedb import MaveDBClient
from genesis_bio_mcp.clients.uniprot import _parse_entry
from genesis_bio_mcp.models import SimilarCompounds


@respx.mock
async def test_mavedb_filters_to_queried_gene() -> None:
    payload = {
        "scoreSets": [
            {
                "urn": "urn:mavedb:00000001-a-1",
                "title": "BRCA1 SGE",
                "numVariants": 100,
                "targetGenes": [{"name": "BRCA1", "uniprotIdFromMappedMetadata": "P38398"}],
            },
            {
                "urn": "urn:mavedb:00000002-a-1",
                "title": "BAP1 SGE (abstract mentions BRCA1)",
                "numVariants": 200,
                "targetGenes": [{"name": "BAP1", "uniprotIdFromMappedMetadata": "Q96TC6"}],
            },
        ]
    }
    respx.post(url__regex=r".*score-sets/search.*").mock(
        return_value=httpx.Response(200, json=payload)
    )
    async with httpx.AsyncClient() as http:
        result = await MaveDBClient(http)._fetch("BRCA1")
    assert result is not None
    assert result.total_score_sets == 1  # the BAP1 set is filtered out
    assert result.score_sets[0].target_gene == "BRCA1"


def test_uniprot_variant_parsing_skips_empty_and_extracts_disease() -> None:
    entry = {
        "primaryAccession": "P04637",
        "entryType": "UniProtKB reviewed (Swiss-Prot)",
        "features": [
            {
                "type": "Natural variant",
                "location": {"start": {"value": 175}},
                "alternativeSequence": {"originalSequence": "R", "alternativeSequences": ["H"]},
                "description": "In LFS; somatic mutation; loss of function.",
            },
            {
                "type": "Natural variant",
                "location": {"start": {"value": 521}},
                "alternativeSequence": {"originalSequence": "R", "alternativeSequences": ["K"]},
                "description": "",
            },
            {
                "type": "Natural variant",
                "location": {"start": {"value": 30}},
                "alternativeSequence": {},  # no substitution -> skipped (was a bare "30")
            },
        ],
    }
    info = _parse_entry(entry, "TP53")
    by_pos = {v.position: v for v in info.known_variants}
    assert "30" not in by_pos  # empty feature skipped
    assert by_pos["175"].disease == "LFS"  # disease token extracted
    assert by_pos["521"].disease is None  # polymorphism: no invented disease
    # And it never renders the misleading "unknown disease" label.
    assert "unknown disease" not in info.to_markdown()


def test_similar_compounds_header_distinguishes_retrieval_from_ranking_metric() -> None:
    md = SimilarCompounds(
        query_smiles="CC(=O)Oc1ccccc1C(=O)O",
        mode="similarity",
        threshold=0.85,
        total_found=0,
        hits=[],
    ).to_markdown()
    assert "Tanimoto ≥ 0.85" not in md  # the false claim is gone
    assert "RDKit ECFP4" in md  # the displayed metric is named honestly
