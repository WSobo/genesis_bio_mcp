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

- Register with `@mcp.tool(annotations=ToolAnnotations(...))` — all four hint fields required
- Input model: inherit `_GeneInput` for single-gene tools; use full `ConfigDict` for multi-field inputs
- Always resolve aliases: `symbol, _ = await _resolve_symbol(params.gene_symbol)`
- Output via `_fmt(result, params.response_format, error_msg)` — never format manually
- Write the docstring for an AI agent: explain *when* to call this tool and *what* the output contains
- Add the client to `server.state` in `lifespan()` and update the tool count in the module docstring

```python
class GetMyDataInput(_GeneInput):
    """Input for get_my_data."""


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
    )
)
async def get_my_data(params: GetMyDataInput) -> str:
    """Fetch data from MyApi for a gene symbol.

    Use this when you need X. Returns markdown with Y and Z fields.

    Args:
        params (GetMyDataInput): gene_symbol, response_format.
    """
    symbol, _ = await _resolve_symbol(params.gene_symbol)
    result = await mcp.state.my_api.get_data(symbol)
    return _fmt(result, params.response_format, f"No MyApi data found for '{symbol}'.")
```

### 4. Tests (`tests/test_clients.py`)

- Use `respx` to mock HTTP at the transport layer — never make real network calls in unit tests
- Test the happy path, a 404/empty response, and an error response
- Verify the returned model has the expected fields
- Env var mocking: `monkeypatch.setenv` / `monkeypatch.delenv(raising=False)` — never mutate `os.environ`
- Add `_MOCK_*` response dicts as module-level constants at the top of the test block

```python
_MOCK_MY_API_RESPONSE = {"gene": "BRAF", "value": 42}


@respx.mock
async def test_my_api_get_data(http_client):
    respx.get(url__regex=r"myapi\.example\.com/data").mock(
        return_value=httpx.Response(200, json=_MOCK_MY_API_RESPONSE)
    )
    client = MyApiClient(http_client)
    result = await client.get_data("BRAF")
    assert result is not None
    assert result.gene_symbol == "BRAF"
```

### 5. Settings (`src/genesis_bio_mcp/config/settings.py`)

If the new client has a rate-limit semaphore or tunable timeout, add a field to the `Settings` class:

```python
my_api_semaphore_limit: int = Field(
    default=3,
    description="Max concurrent requests to MyApi.",
    gt=0,
)
```

Then in the client:
```python
from genesis_bio_mcp.config.settings import settings
_SEMAPHORE = asyncio.Semaphore(settings.my_api_semaphore_limit)
```

Users can override via `GENESIS_MY_API_SEMAPHORE_LIMIT=5` or a `.env` file.

### 6. Workflow agent (`src/genesis_bio_mcp/workflow_agent.py`)

Every new tool must also be reachable from `run_biology_workflow`:

- Add `async def _get_<name>_fn(...)` inside `build_tool_registry()`
- Add a `ToolSpec(...)` entry: `tool_category`, `use_when` (one embedding-searchable sentence), `input_schema`, `fn`
- In `tests/test_workflow_agent.py`: add client attr name to `_mock_state()`, add client method to the inner loop, add tool name to `expected_tools`, bump the `"N tools"` assertion count

### 7. Checklist before opening a PR

- [ ] `uv run ruff format . && uv run ruff check --fix .` — no errors
- [ ] `uv run pytest tests/ -v` — all passing
- [ ] README tools table updated
- [ ] API rate limit and auth notes added to README API reference table
- [ ] If new client has semaphore/timeout: added field to `config/settings.py`
- [ ] `CHANGELOG.md` updated under `## [Unreleased]`
- [ ] `workflow_agent.py` tool registry updated and `test_workflow_agent.py` passing

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
