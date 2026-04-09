# genesis-bio-mcp

An MCP (Model Context Protocol) server that connects AI agents to major public biomedical databases for drug discovery target prioritization.

Ask an AI agent *"Assess BRAF as an oncology target for melanoma"* and watch it autonomously chain queries across UniProt, Open Targets, DepMap, GWAS Catalog, and PubChem into a structured evidence report — no hardcoded scripts, no manual API calls.

## Tools

| Tool | Database | Purpose |
|------|----------|---------|
| `resolve_gene` | UniProt + NCBI | Resolve gene aliases → canonical HGNC symbol, NCBI ID, UniProt accession |
| `get_protein_info` | UniProt Swiss-Prot | Protein function, pathways, disease variants, PDB structures |
| `get_target_disease_association` | Open Targets | Evidence-based association score (0–1) for a target–disease pair |
| `get_cancer_dependency` | DepMap | CRISPR essentiality scores across cancer cell lines |
| `get_gwas_evidence` | GWAS Catalog | Genome-wide significant SNP associations for a trait |
| `get_compounds` | PubChem | Active small molecules with potency data (IC50/EC50) |
| `prioritize_target` | All of the above | Full parallel evidence synthesis → priority score (0–10) + report |

## Quickstart

### Install

```bash
pip install genesis-bio-mcp
```

Or from source:

```bash
git clone https://github.com/WSobo/genesis-bio-mcp
cd genesis-bio-mcp
pip install -e ".[dev]"
```

### Add to Claude Desktop

Edit `~/Library/Application Support/Claude/claude_desktop_config.json` (macOS) or `%APPDATA%\Claude\claude_desktop_config.json` (Windows):

```json
{
  "mcpServers": {
    "genesis-bio-mcp": {
      "command": "uvx",
      "args": ["genesis-bio-mcp"],
      "env": {
        "NCBI_EMAIL": "your@email.com"
      }
    }
  }
}
```

Restart Claude Desktop. You'll see the genesis-bio-mcp tools available in your conversation.

### Try it

Ask Claude:
- *"Is PCSK9 a good target for cardiovascular disease? Check Open Targets and PubChem."*
- *"Assess BRAF as an oncology target for melanoma — full report."*
- *"What GWAS evidence links FTO to obesity?"*
- *"How many active compounds exist for EGFR in PubChem?"*

## Example Output

```
prioritize_target("BRAF", "melanoma")

→ priority_score: 8.4 / 10
→ priority_tier: High
→ evidence_summary: "BRAF shows strong Open Targets association with melanoma (score: 0.89,
  n=312 evidence items). DepMap CRISPR data show dependency in 61% of cancer lines, highest
  in Skin, Thyroid, Colon. GWAS Catalog links 5 variants near BRAF to melanoma-related traits
  (strongest p=1.00e-15). PubChem reports 312 active compounds against BRAF, indicating strong
  druggability."

→ disease_association.overall_score: 0.89
→ disease_association.somatic_mutation_score: 0.95
→ disease_association.known_drug_score: 0.88
→ cancer_dependency.fraction_dependent_lines: 0.61
→ compounds.total_active_compounds: 312
→ compounds[0]: Dabrafenib, IC50 = 6.0 nM
```

## Architecture

```
src/genesis_bio_mcp/
├── server.py                  # FastMCP server, tool registration, shared httpx client
├── models.py                  # Pydantic output models (all fields documented for agents)
├── clients/
│   ├── uniprot.py             # UniProt REST: gene_exact query, Swiss-Prot parsing
│   ├── open_targets.py        # Open Targets GraphQL: 3-step resolution (gene→Ensembl, disease→EFO, assoc)
│   ├── depmap.py              # DepMap portal API with graceful fallback
│   ├── gwas.py                # GWAS Catalog HAL/REST, Unicode normalization, trait filtering
│   └── pubchem.py             # PubChem REST, asyncio.Semaphore rate limiting, tenacity retries
└── tools/
    ├── gene_resolver.py       # Multi-source alias resolution (UniProt primary, NCBI E-utils for gene ID)
    └── target_prioritization.py  # asyncio.gather orchestration, safe_call error isolation, score computation
```

**Key design decisions:**
- Single shared `httpx.AsyncClient` via FastMCP `lifespan` for connection pooling
- `asyncio.gather` in `prioritize_target` for parallel database queries (5 APIs simultaneously)
- Every sub-query wrapped in `safe_call` — agent never crashes on a single API failure
- Opinionated output filtering — agents see 8 meaningful fields, not raw 50-field API blobs
- Agent-readable tool docstrings that explain *when* to use each tool and *how* to format inputs

## Development

```bash
pip install -e ".[dev]"
pytest tests/ -v
pytest tests/ -v --cov=genesis_bio_mcp
```

### Running the server directly

```bash
genesis-bio-mcp          # stdio transport (for MCP clients)
python -m genesis_bio_mcp.server
```

### Environment variables

| Variable | Default | Description |
|----------|---------|-------------|
| `NCBI_EMAIL` | `genesis-bio-mcp@example.com` | Required for NCBI E-utils polite use |

## API Notes

| Database | API Type | Rate Limit | Notes |
|----------|----------|------------|-------|
| UniProt | REST | Generous | Filter `organism_id:9606` + `reviewed:true` for human Swiss-Prot |
| Open Targets | GraphQL | Generous (slow ~2s) | Requires Ensembl ID and EFO ID — resolved automatically |
| DepMap | Undocumented REST | N/A | No stable public API; graceful fallback with portal link |
| GWAS Catalog | REST/HAL | Moderate | HAL JSON; Unicode normalization required for trait matching |
| PubChem | REST | 5 req/sec | Returns HTTP 503 on rate limit; handled with tenacity + Semaphore |
| NCBI E-utils | REST | 3 req/sec | Requires `email` parameter; set `NCBI_EMAIL` env var |

## License

MIT
