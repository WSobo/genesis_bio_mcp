# Genesis Bio MCP Rules

## Package Management
PKG: `uv` ONLY. NO pip/conda/venv. Add: `uv add [--dev] <pkg>`. Run: `uv run <cmd>`.

## Filesystem & Search
NO built-in Read/Grep/Glob. NO raw ls/cat/git. MUST use bash `rtk <cmd>` (`rtk ls`, `rtk read <file>`, `rtk grep`, `rtk git`).

## Linting
`uv run ruff format .` & `uv run ruff check --fix .` (ALWAYS run before commit).

## Testing
`uv run pytest tests/ -v`

## Git
- Use `rtk git` for ALL version control.
- Commits MUST be atomic (one feature/fix per commit).
- Conventional Commits ONLY (`feat:`, `fix:`, `refactor:`, `docs:`).
- NEVER `git push` without explicit user permission.
- NEVER track/commit `data/`, `.parquet`, or `.csv` files.

---

## Architecture

### Invariants
- Type hints MANDATORY everywhere.
- Pydantic V2 ONLY (`model_dump`, `model_validate`, `model_dump_json`).
- MCP tools output strictly formatted MARKDOWN strings, NEVER raw JSON dicts (unless `response_format="json"`).
- All tools must include `ToolAnnotations` (`readOnlyHint`, `destructiveHint`, `idempotentHint`, `openWorldHint`).
- Tool names: `{service}_{action}_{resource}` snake_case.
- All tools support `response_format` param (`"markdown"` default, `"json"` for programmatic use).

### Client pattern (`src/genesis_bio_mcp/clients/<source>.py`)

- `__init__` accepts a single `httpx.AsyncClient` stored as `self._client` — never create clients internally.
- Module-level semaphore: `_SEMAPHORE = asyncio.Semaphore(N)` — N=2 for aggressive APIs (BioGRID, STRING), N=3 for most others. Acquire with `async with _SEMAPHORE:`.
- Optional session cache: `self._cache: dict[str, Model] = {}` — check before calling, store only on success. NEVER cache failures or `None` results.
- All HTTP calls wrapped in `try/except Exception` — return `None` or `[]` on error; `logger.warning("API failed for %s: %s", symbol, exc)`. Never raise to the caller.
- Always set explicit `timeout=` on every HTTP call (20–25 s is standard).
- Key-guard: check `os.environ.get("API_KEY")` at the top of the public method; return `None` immediately with a `logger.warning` if absent.
- `logger.debug(...)` for cache hits; `logger.warning(...)` for recoverable failures.

```python
_SEMAPHORE = asyncio.Semaphore(3)

class MyApiClient:
    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client
        self._cache: dict[str, MyModel] = {}

    async def get_data(self, gene_symbol: str) -> MyModel | None:
        symbol = gene_symbol.strip().upper()
        if symbol in self._cache:
            logger.debug("MyApi cache hit: %s", symbol)
            return self._cache[symbol]
        async with _SEMAPHORE:
            try:
                resp = await self._client.get(URL, timeout=20.0)
                resp.raise_for_status()
                result = _parse(resp.json(), symbol)
            except Exception as exc:
                logger.warning("MyApi failed for %s: %s", symbol, exc)
                return None
        self._cache[symbol] = result
        return result
```

### Model pattern (`src/genesis_bio_mcp/models.py`)

- Inherit from `pydantic.BaseModel`.
- Every model **must** implement `to_markdown(self) -> str` — human-readable, agent-parseable.
- All fields: `Field(description="...")` — required for MCP tool schema auto-generation.
- Optional fields: `foo: str | None = Field(None, description="...")`.
- Group models by data source; add to the appropriate section in `models.py`.
- Pydantic V2 only — no V1 patterns (`parse_obj`, `schema()`, validators as class methods without decorators).

### Tool registration pattern (`src/genesis_bio_mcp/server.py`)

1. **Input model** — inherit `_GeneInput` for single-gene tools (provides `gene_symbol` + `response_format`):
   ```python
   class GetMyDataInput(_GeneInput):
       """Input for get_my_data."""
   ```
   Multi-field inputs: `ConfigDict(str_strip_whitespace=True, validate_assignment=True, extra="forbid")`.

2. **Tool function**:
   ```python
   @mcp.tool(
       annotations=ToolAnnotations(
           readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True
       )
   )
   async def get_my_data(params: GetMyDataInput) -> str:
       """Describe when to call this and what the output contains — written for an AI agent."""
       symbol, _ = await _resolve_symbol(params.gene_symbol)
       result = await mcp.state.my_api.get_data(symbol)
       return _fmt(result, params.response_format, f"No data found for '{symbol}'.")
   ```

3. **Wire the client** in `lifespan()`: `server.state.my_api = MyApiClient(client)`.
4. **Update the tool count** in the module docstring at the top of `server.py`.

Key rules:
- Always resolve aliases via `_resolve_symbol()` before querying any database.
- Always use `_fmt(result, params.response_format, error_msg)` — never format output manually.
- Docstrings are the tool's API contract with the LLM — explain *when* and *what*, not *how*.

### Workflow agent pattern (`src/genesis_bio_mcp/workflow_agent.py`)

Every new MCP tool must also be available inside `run_biology_workflow`:

1. Add `async def _get_<name>_fn(...)` inside `build_tool_registry()`.
2. Add a `ToolSpec(...)` entry to the returned dict:
   - `tool_category`: one of `gene_annotation`, `disease_evidence`, `druggability`, `structure`, `pathways`, `synthesis`
   - `use_when`: a single embedding-searchable sentence describing when to invoke the tool
   - `input_schema`: JSON Schema dict matching the tool's required inputs
   - `fn`: the inner function
3. In `tests/test_workflow_agent.py`:
   - Add client attr to `_mock_state()` attrs list and add its method to the inner method loop
   - Add tool name to `expected_tools` in `test_build_tool_registry_has_all_tools`
   - Bump `"N tools"` assertion in `test_format_registry_docs_structure`

### Testing pattern (`tests/`)

- HTTP mocking: `@respx.mock` decorator; `respx.get(url__regex=r"...")` for dynamic URLs.
- Shared fixtures and widely-reused mock data live in `tests/conftest.py` as `MOCK_*` constants.
- Client-local mock data: `_MOCK_*` module-level dicts at the top of `tests/test_clients.py`.
- Required coverage per client: **happy path**, **404/empty response**, **error/network failure**.
- Env var mocking: `monkeypatch.setenv("KEY", "val")` / `monkeypatch.delenv("KEY", raising=False)`.
- Never make live network calls in unit tests.

---

## References

- `CONTRIBUTING.md` — step-by-step walkthrough for adding a new data source (5 steps + checklist)
- `src/genesis_bio_mcp/clients/biogrid.py` — canonical new-client example (key-guard, semaphore, cache)
- `src/genesis_bio_mcp/server.py` — tool registration reference (ToolAnnotations, _GeneInput, _fmt, lifespan)
- `src/genesis_bio_mcp/workflow_agent.py` — ToolSpec pattern and build_tool_registry() reference
