# Deployment Guide

## Environment variables

All runtime settings use the `GENESIS_` prefix. Every variable has a sensible default; you only need to set the ones you want to change.

| Variable | Default | Description |
|----------|---------|-------------|
| `ANTHROPIC_API_KEY` | *(required for `run_biology_workflow`)* | Anthropic API key — must be in the MCP server's env, not just the shell |
| `GENESIS_HTTPX_TIMEOUT` | `30.0` | Default HTTP timeout in seconds for all outbound API calls |
| `GENESIS_DEPMAP_CACHE_PATH` | `data/depmap_cache.csv` | Path to DepMap disk cache (multi-MB CSV, refreshed every 7 days) |
| `GENESIS_DEPMAP_CACHE_MAX_AGE_DAYS` | `7` | Days before DepMap cache is re-downloaded |
| `GENESIS_DEPMAP_TASK_TIMEOUT_SECS` | `120.0` | Timeout when polling DepMap async Celery tasks |
| `GENESIS_GWAS_CACHE_PATH` | `data/gwas_cache.json` | Path to GWAS Catalog result cache |
| `GENESIS_GWAS_CACHE_TTL_SECS` | `86400` | GWAS cache TTL in seconds (default 24 h) |
| `GENESIS_EFO_CACHE_PATH` | `data/efo_cache.json` | Path to EFO ontology resolution cache |
| `GENESIS_EFO_CACHE_TTL_SECS` | `604800` | EFO cache TTL in seconds (default 7 days) |
| `GENESIS_CHEMBL_SEMAPHORE_LIMIT` | `2` | Max concurrent requests to ChEMBL |
| `GENESIS_PUBCHEM_SEMAPHORE_LIMIT` | `3` | Max concurrent requests to PubChem |
| `GENESIS_REACTOME_SEMAPHORE_LIMIT` | `3` | Max concurrent requests to Reactome |
| `GENESIS_CLAUDE_MODEL` | `claude-sonnet-4-6` | Claude model used by `run_biology_workflow` |
| `NCBI_EMAIL` | `genesis-bio-mcp@example.com` | Email sent in NCBI E-utils requests (required by NCBI ToS for production) |
| `BIOGRID_ACCESS_KEY` | *(unset)* | Required for `get_biogrid_interactions`. Free key from <https://webservice.thebiogrid.org/>. When unset, the tool returns an empty result with a warning logged — all other tools work normally. |
| `OPENFDA_API_KEY` | *(unset)* | Optional — lifts the OpenFDA free-tier quota (240 req/min, 1000 req/day). Free key from <https://open.fda.gov/apis/authentication/>. Affects `get_drug_history` and `prioritize_target` (extended mode); without it the tools work but may rate-limit on heavy workflows. |
| `GENESIS_OPENFDA_CACHE_PATH` | `data/openfda_cache.json` | Path to OpenFDA drug-safety cache |
| `GENESIS_OPENFDA_CACHE_TTL_SECS` | `604800` | OpenFDA cache TTL in seconds (default 7 days) |

Set via a `.env` file in the working directory (loaded automatically) or as shell exports.

---

## Local / development

```bash
uv sync
export ANTHROPIC_API_KEY=sk-ant-...
uv run genesis-bio-mcp
```

Or with custom settings:
```bash
GENESIS_HTTPX_TIMEOUT=60 GENESIS_CLAUDE_MODEL=claude-opus-4-6 uv run genesis-bio-mcp
```

---

## Claude Desktop

Edit `~/Library/Application Support/Claude/claude_desktop_config.json` (macOS) or `%APPDATA%\Claude\claude_desktop_config.json` (Windows):

```json
{
  "mcpServers": {
    "genesis-bio-mcp": {
      "command": "uv",
      "args": ["run", "--directory", "/absolute/path/to/genesis-bio-mcp", "genesis-bio-mcp"],
      "env": {
        "ANTHROPIC_API_KEY": "sk-ant-...",
        "GENESIS_DEPMAP_CACHE_PATH": "/absolute/path/to/genesis-bio-mcp/data/depmap_cache.csv",
        "GENESIS_GWAS_CACHE_PATH": "/absolute/path/to/genesis-bio-mcp/data/gwas_cache.json",
        "GENESIS_EFO_CACHE_PATH": "/absolute/path/to/genesis-bio-mcp/data/efo_cache.json"
      }
    }
  }
}
```

> **Note:** `ANTHROPIC_API_KEY` must be in the `env` block here — the shell environment where Claude Desktop launches is not forwarded to MCP server processes. If `run_biology_workflow` returns an auth error, this is the most common cause.

Use absolute paths for cache files when running via Claude Desktop, since the working directory is not predictable.

---

## Docker

```dockerfile
FROM python:3.11-slim

# Install uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

WORKDIR /app
COPY pyproject.toml uv.lock ./
COPY src/ ./src/

RUN uv sync --no-dev

# Cache files live in a volume so they persist across container restarts
VOLUME ["/app/data"]

ENV GENESIS_DEPMAP_CACHE_PATH=/app/data/depmap_cache.csv \
    GENESIS_GWAS_CACHE_PATH=/app/data/gwas_cache.json \
    GENESIS_EFO_CACHE_PATH=/app/data/efo_cache.json

# ANTHROPIC_API_KEY must be injected at runtime — do not bake into image
CMD ["uv", "run", "genesis-bio-mcp"]
```

```bash
docker build -t genesis-bio-mcp .
docker run \
  -e ANTHROPIC_API_KEY=sk-ant-... \
  -v $(pwd)/data:/app/data \
  -p 3000:3000 \
  genesis-bio-mcp
```

The `data/` volume is important: without it, DepMap's ~10 MB CSV and the GWAS/EFO caches are re-downloaded on every container start.

---

## Production checklist

- [ ] Set `ANTHROPIC_API_KEY` in the server environment (not just the host shell)
- [ ] Set `NCBI_EMAIL` to a real address — NCBI rate-limits anonymous requests
- [ ] Mount a persistent volume for `data/` cache files
- [ ] Use absolute paths for all `GENESIS_*_CACHE_PATH` variables
- [ ] Consider increasing `GENESIS_HTTPX_TIMEOUT` (to 60s) on slow networks
- [ ] Pin `GENESIS_CLAUDE_MODEL` to a specific model ID for reproducibility
