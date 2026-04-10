"""Full integration test — runs prioritize_target against live APIs for all test cases.

Usage:
    uv run python test_full.py
    uv run python test_full.py --gene BRAF --disease melanoma   # single run

Outputs:
    - Console: Markdown report + score summary for each case
    - examples/<gene>_<disease>_report.json  (raw JSON for debugging)
    - examples/<gene>_<disease>_report.md    (Markdown for README reference)
"""

from __future__ import annotations

import asyncio
import json
import re
import sys
import time
from pathlib import Path

import httpx

from genesis_bio_mcp.clients.chembl import ChEMBLClient
from genesis_bio_mcp.clients.depmap import DepMapClient, load_depmap_cache
from genesis_bio_mcp.clients.gwas import GwasClient
from genesis_bio_mcp.clients.open_targets import OpenTargetsClient
from genesis_bio_mcp.clients.pubchem import PubChemClient
from genesis_bio_mcp.clients.uniprot import UniProtClient
from genesis_bio_mcp.tools.target_prioritization import prioritize_target

HEADERS = {
    "User-Agent": "genesis-bio-mcp/0.1 (research; github.com/WSobo/genesis-bio-mcp)",
    "Accept": "application/json",
}

# ---------------------------------------------------------------------------
# Test matrix
# ---------------------------------------------------------------------------

TEST_CASES: list[tuple[str, str]] = [
    # (gene_input, indication)              # what this tests
    ("BRAF", "melanoma"),  # baseline — classic oncogene, activating V600E
    (
        "EGFR",
        "non-small cell lung carcinoma",
    ),  # targeted therapy archetype, massive compound data
    ("KRAS", "pancreatic cancer"),  # historically "undruggable", high DepMap dependency
    ("HER2", "breast cancer"),  # resolve_gene: "HER2" → ERBB2
    ("PCSK9", "hypercholesterolemia"),  # GWAS success story, biologic drugs
    ("FTO", "obesity"),  # massive GWAS hits, complex trait genetics
    ("TNF", "rheumatoid arthritis"),  # autoimmune, heavy biologic/compound data
    ("PTGS2", "inflammation"),  # COX-2: flood of NSAIDs/Coxibs in PubChem
    ("TP53", "squamous cell carcinoma"),  # most-mutated cancer gene, low druggability
    ("CD274", "melanoma"),  # PD-L1 checkpoint (alias resolution)
    ("p53", "lung cancer"),  # synonym → TP53
    ("COX2", "pain"),  # alias → PTGS2
]


def _safe_filename(s: str) -> str:
    return re.sub(r"[^a-zA-Z0-9_-]", "_", s).lower()


async def run_one(
    gene: str,
    indication: str,
    *,
    uniprot: UniProtClient,
    open_targets: OpenTargetsClient,
    depmap: DepMapClient,
    gwas: GwasClient,
    pubchem: PubChemClient,
    chembl: ChEMBLClient,
    examples_dir: Path,
) -> dict:
    t0 = time.perf_counter()
    print(f"\n{'=' * 60}")
    print(f"  {gene!r} × {indication!r}")
    print(f"{'=' * 60}")

    result = await prioritize_target(
        gene_symbol=gene,
        indication=indication,
        uniprot=uniprot,
        open_targets=open_targets,
        depmap=depmap,
        gwas=gwas,
        pubchem=pubchem,
        chembl=chembl,
    )
    elapsed = time.perf_counter() - t0

    resolved = result.resolution.hgnc_symbol
    if resolved != gene.upper():
        print(f"  Resolved: {gene!r} → {resolved}")

    print(f"  Score:    {result.priority_score:.1f}/10  [{result.priority_tier}]")
    print(f"  Gaps:     {result.data_gaps or 'none'}")
    print(f"  Errors:   {list(result.errors.keys()) or 'none'}")

    if result.cancer_dependency:
        src = (
            "DepMap (real)"
            if "DepMap Chronos" in result.cancer_dependency.data_source
            else "OT proxy"
        )
        print(
            f"  DepMap:   {src} — {int(result.cancer_dependency.fraction_dependent_lines * 100)}% dependent"
        )

    print(f"  Time:     {elapsed:.1f}s")
    print()
    print(result.evidence_summary)

    # Save files
    stem = f"{_safe_filename(gene)}_{_safe_filename(indication)}"
    json_path = examples_dir / f"{stem}_report.json"
    md_path = examples_dir / f"{stem}_report.md"

    with open(json_path, "w") as f:
        json.dump(result.model_dump(), f, indent=2, default=str)

    with open(md_path, "w") as f:
        f.write(result.to_markdown())

    return {
        "gene_input": gene,
        "resolved": resolved,
        "indication": indication,
        "score": result.priority_score,
        "tier": result.priority_tier,
        "data_gaps": result.data_gaps,
        "errors": list(result.errors.keys()),
        "depmap_real": result.cancer_dependency is not None
        and "DepMap Chronos" in result.cancer_dependency.data_source,
        "elapsed_s": round(elapsed, 1),
    }


async def run(cases: list[tuple[str, str]]) -> None:
    examples_dir = Path("examples")
    examples_dir.mkdir(exist_ok=True)

    async with httpx.AsyncClient(headers=HEADERS, timeout=60.0, follow_redirects=True) as http:
        print("Loading DepMap gene_dep_summary cache...")
        cache = await load_depmap_cache(http)
        print(
            f"  {len(cache)} genes loaded from DepMap"
            if cache
            else "  Cache empty — will use OT proxy"
        )

        clients = dict(
            uniprot=UniProtClient(http),
            open_targets=OpenTargetsClient(http),
            depmap=DepMapClient(http, cache),
            gwas=GwasClient(http),
            pubchem=PubChemClient(http),
            chembl=ChEMBLClient(http),
        )

        summary: list[dict] = []
        for gene, indication in cases:
            row = await run_one(gene, indication, examples_dir=examples_dir, **clients)
            summary.append(row)

    # Print summary table
    print(f"\n\n{'=' * 70}")
    print("SUMMARY")
    print(f"{'=' * 70}")
    print(
        f"{'Input':<10} {'Resolved':<10} {'Indication':<30} {'Score':>6} {'Tier':<8} {'DepMap':>8} {'t(s)':>5}"
    )
    print("-" * 70)
    for r in summary:
        dep = "real" if r["depmap_real"] else "proxy"
        gaps = ",".join(r["data_gaps"]) if r["data_gaps"] else "-"
        print(
            f"{r['gene_input']:<10} {r['resolved']:<10} {r['indication'][:29]:<30} "
            f"{r['score']:>5.1f}  {r['tier']:<8} {dep:>8} {r['elapsed_s']:>4.1f}s"
            + (f"  [gaps: {gaps}]" if r["data_gaps"] else "")
        )

    json_summary = Path("examples") / "test_summary.json"
    with open(json_summary, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSummary saved to {json_summary}")


if __name__ == "__main__":
    # Allow single-case override: python test_full.py BRAF melanoma
    args = sys.argv[1:]
    if len(args) == 2:
        cases = [(args[0], args[1])]
    elif "--gene" in args and "--disease" in args:
        g = args[args.index("--gene") + 1]
        d = args[args.index("--disease") + 1]
        cases = [(g, d)]
    else:
        cases = TEST_CASES

    asyncio.run(run(cases))
