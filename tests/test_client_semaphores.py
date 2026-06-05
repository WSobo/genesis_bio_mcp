"""Invariant: every direct-HTTP client module exposes a concurrency semaphore.

These six clients previously issued unbounded concurrent HTTP fan-out (no
module-level ``_SEMAPHORE``), unlike the documented client pattern and the
other 18 clients. This locks the fix in so the pattern can't silently regress.
"""

from __future__ import annotations

import asyncio
import importlib

import pytest

# Direct-HTTP clients that were missing a module semaphore before this fix.
_MODULES = [
    "alphafold",
    "clinical_trials",
    "dgidb",
    "gwas",
    "open_targets",
    "uniprot",
]


@pytest.mark.parametrize("name", _MODULES)
def test_client_module_has_semaphore(name: str) -> None:
    mod = importlib.import_module(f"genesis_bio_mcp.clients.{name}")
    sem = getattr(mod, "_SEMAPHORE", None)
    assert isinstance(sem, asyncio.Semaphore), f"{name} is missing a module-level _SEMAPHORE"
