## What

<!-- One-paragraph summary of what this PR does. -->

## Why

<!-- Motivation: what problem does this solve, or what feature does it add? -->

## How

<!-- Brief description of the approach: key design decisions, trade-offs. -->

## Checklist

- [ ] New client added to `clients/` with `_safe()` wrapping
- [ ] New Pydantic V2 models added to `models.py` with `to_markdown()`
- [ ] Tool registered in `server.py` with agent-readable docstring
- [ ] Unit tests added (mock at HTTP layer with `respx`)
- [ ] `uv run ruff format . && uv run ruff check --fix .` passed
- [ ] `uv run pytest tests/ -v` passes (all green)
- [ ] README updated if new tool or scoring change

## API Endpoints Touched

<!-- List any external APIs added or modified. Include rate limit and auth notes. -->

| API | URL | Auth | Rate Limit |
|-----|-----|------|-----------|
|     |     |      |           |
