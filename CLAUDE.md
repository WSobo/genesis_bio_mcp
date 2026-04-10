# Genesis Bio MCP - Development Rules

## Package Management
- This project strictly uses `uv`. Do NOT use `pip`, `conda`, or `venv`.
- Add dependencies: `uv add <package>` (use `--dev` for dev dependencies).
- Execute scripts: `uv run <command>`.

## Token Efficiency & File Reading
- **CRITICAL:** Do NOT use your built-in `Read`, `Grep`, or `Glob` tools.
- To read files, view directories, or search code, you MUST explicitly use the `rtk` CLI utility via the Bash tool.
- NEVER use raw `ls`, `cat`, or `git` commands.
- Use `rtk ls` instead of `ls`.
- Use `rtk read <file>` instead of `cat`.
- Use `rtk grep <pattern>` instead of `grep`.
- Use `rtk git <cmd>` instead of `git <cmd>`.

## Architecture & Coding Style
- We use **Pydantic V2**. Strictly use V2 syntax (`model_dump()` instead of `dict()`, `model_validate()` instead of `parse_obj()`).
- All MCP tools must return strictly formatted **Markdown strings** (not raw JSON dicts) to ensure downstream LLMs can easily parse the reports.
- Wrap all external API calls in `safe_call` blocks to prevent the whole parallel pipeline from failing if one database timeouts.

## Testing
- Run the test suite: `uv run pytest tests/ -v`

