import asyncio
import httpx
from genesis_bio_mcp.clients.uniprot import UniProtClient
from genesis_bio_mcp.clients.open_targets import OpenTargetsClient
from genesis_bio_mcp.clients.depmap import DepMapClient
from genesis_bio_mcp.clients.gwas import GwasClient
from genesis_bio_mcp.clients.pubchem import PubChemClient
from genesis_bio_mcp.tools.target_prioritization import prioritize_target
import json

async def run():
    async with httpx.AsyncClient(
        headers={"User-Agent": "genesis-bio-mcp/0.1 (research; github.com/WSobo/genesis-bio-mcp)"},
        timeout=60.0
    ) as http:
        result = await prioritize_target(
            gene_symbol="BRAF",
            indication="melanoma",
            uniprot=UniProtClient(http),
            open_targets=OpenTargetsClient(http),
            depmap=DepMapClient(http),
            gwas=GwasClient(http),
            pubchem=PubChemClient(http),
        )
        print(f"Score:   {result.priority_score}")
        print(f"Tier:    {result.priority_tier}")
        print(f"Gaps:    {result.data_gaps}")
        print(f"\nSummary:\n{result.evidence_summary}")
        
        # Save full JSON output for README
        with open("examples/braf_melanoma_report.json", "w") as f:
            json.dump(result.model_dump(), f, indent=2, default=str)
        print("\nFull report saved to examples/braf_melanoma_report.json")

asyncio.run(run())
