# Contributing to genesis-bio-mcp

## Setup

```bash
git clone https://github.com/WSobo/genesis-bio-mcp
cd genesis-bio-mcp
uv sync
```

## Development workflow

```bash
# Run tests
uv run pytest tests/ -v

# Lint and format (always before committing)
uv run ruff format .
uv run ruff check --fix .
```

## Adding a new data source

Follow this pattern for every new API integration:

### 1. Client (`src/genesis_bio_mcp/clients/<source>.py`)

- Accept a single `httpx.AsyncClient` in `__init__`
- Wrap all HTTP calls in `try/except Exception` — return `None` or `[]` on failure, log with `logger.warning`
- Apply `asyncio.Semaphore` if the API has rate limits (see table in README)
- Never raise exceptions to the caller

```python
class MyApiClient:
    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client

    async def get_data(self, gene_symbol: str) -> Optional[MyModel]:
        try:
            resp = await self._client.get(URL, timeout=20.0)
            resp.raise_for_status()
            return _parse(resp.json())
        except Exception as exc:
            logger.warning("MyApi failed for %s: %s", gene_symbol, exc)
            return None
```

### 2. Model (`src/genesis_bio_mcp/models.py`)

- Use Pydantic V2 (`BaseModel`, `model_validator`, `Field`)
- Every model needs a `to_markdown(self) -> str` method
- Output must be human-readable and agent-parseable — avoid deeply nested structures

```python
class MyModel(BaseModel):
    gene_symbol: str
    some_field: Optional[float] = None

    def to_markdown(self) -> str:
        lines = [f"## MyData — {self.gene_symbol}"]
        if self.some_field is not None:
            lines.append(f"- Some field: {self.some_field:.2f}")
        return "\n".join(lines)
```

### 3. Tool (`src/genesis_bio_mcp/server.py`)

- Register with `@mcp.tool()`
- The function **must return a markdown string** — never a dict, never raw JSON
- Write the docstring for an AI agent: explain *when* to call this tool and *what* the inputs mean
- Inject the client from `server.state`

```python
@mcp.tool()
async def get_my_data(gene_symbol: str) -> str:
    """Fetch data from MyApi for a gene symbol.

    Use this when you need X. Input must be a canonical HGNC gene symbol (e.g. BRAF, EGFR).
    Returns markdown with Y and Z fields.
    """
    result = await request.app.state.my_api.get_data(gene_symbol)
    if result is None:
        return f"## MyData — {gene_symbol}\n\nNo data found."
    return result.to_markdown()
```

- Add the client to the `lifespan` context manager in `server.py`

### 4. Tests (`tests/test_clients.py`)

- Use `respx` to mock HTTP at the transport layer — never make real network calls in unit tests
- Test the happy path, a 404/empty response, and an error response
- Verify the returned model has the expected fields

```python
@respx.mock
async def test_my_api_get_data(http_client):
    respx.get("https://myapi.example.com/data").mock(
        return_value=httpx.Response(200, json={"result": "value"})
    )
    client = MyApiClient(http_client)
    result = await client.get_data("BRAF")
    assert result is not None
    assert result.gene_symbol == "BRAF"
```

### 5. Checklist before opening a PR

- [ ] `uv run ruff format . && uv run ruff check --fix .` — no errors
- [ ] `uv run pytest tests/ -v` — all passing
- [ ] README tools table updated
- [ ] API rate limit and auth notes added to README API reference table

## Commit style

Conventional Commits only:

```
feat: add Reactome pathway client
fix: handle 204 empty response from RCSB search
refactor: extract _parse_interactions helper in dgidb client
docs: update README with extended mode example
```

One feature or fix per commit. Never `git push` without explicit review.

## MCP tool design principles

- **Markdown output only** — tools are consumed by language models; structured text beats JSON
- **Fail gracefully** — a broken API returns a "no data" message, not an exception traceback
- **Agent-readable docstrings** — the docstring IS the tool's API contract with the LLM
- **Opinionated filtering** — expose 5–10 meaningful fields, not raw API blobs
- **Score stability** — new data sources extend context; they don't silently change scoring without a deliberate versioned change
