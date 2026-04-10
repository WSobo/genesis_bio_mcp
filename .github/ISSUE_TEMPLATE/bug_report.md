---
name: Bug report
about: API failure, wrong output, or crash
labels: bug
---

## Description

<!-- What went wrong? -->

## Tool / Gene / Indication

```
Tool: prioritize_target / get_drug_history / ...
Gene: BRAF
Indication: melanoma
```

## Observed vs Expected

**Observed:**
```
...
```

**Expected:**
```
...
```

## Reproduction

```bash
uv run python test_full.py BRAF melanoma
# or paste the MCP call
```

## Environment

- OS:
- Python version (`python --version`):
- Package version (`uv pip show genesis-bio-mcp`):

## Logs

<!-- Paste relevant log output (set `PYTHONLOGLEVEL=DEBUG` for verbose output) -->
