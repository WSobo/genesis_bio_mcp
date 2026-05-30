"""Server health — upstream reachability + runtime/cache status.

Backs the ``health://status`` MCP resource: lets a client, orchestrator, or ops
process verify the server and its critical upstream services are reachable
before relying on them. Pure helpers here; the resource wiring lives in server.py.
"""

from __future__ import annotations

import asyncio
import time

import httpx

from genesis_bio_mcp import __version__
from genesis_bio_mcp.config.settings import settings


def _upstream_checks() -> list[tuple[str, str]]:
    """(label, cheap-reachability-URL) for a representative set of critical upstreams."""
    return [
        ("UniProt", "https://rest.uniprot.org/uniprotkb/search?query=BRAF&size=1&format=list"),
        ("Open Targets", "https://api.platform.opentargets.org/api/v4/graphql"),
        ("AlphaFold", "https://alphafold.ebi.ac.uk/api/prediction/P15056"),
        (
            "PubChem",
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/"
            "MolecularFormula/JSON",
        ),
        ("Reactome", "https://reactome.org/ContentService/data/database/version"),
        ("STRING", "https://string-db.org/api/json/version"),
        ("UMA-Inverse", f"{settings.uma_api_url.rstrip('/')}/health"),
    ]


async def _probe(client: httpx.AsyncClient, label: str, url: str) -> dict:
    """Probe one upstream; any HTTP response = reachable, a 5xx = degraded."""
    t0 = time.perf_counter()
    try:
        resp = await client.get(url, timeout=8.0)
        ms = round((time.perf_counter() - t0) * 1000)
        return {
            "name": label,
            "reachable": True,
            "status": resp.status_code,
            "latency_ms": ms,
            "ok": resp.status_code < 500,
        }
    except Exception as exc:
        ms = round((time.perf_counter() - t0) * 1000)
        return {
            "name": label,
            "reachable": False,
            "status": None,
            "latency_ms": ms,
            "ok": False,
            "error": type(exc).__name__,
        }


async def check_upstreams(client: httpx.AsyncClient) -> list[dict]:
    """Probe all critical upstreams concurrently and return per-service status."""
    checks = _upstream_checks()
    return list(await asyncio.gather(*[_probe(client, label, url) for label, url in checks]))


def build_health_report(
    *,
    n_tools: int,
    uptime_s: float,
    depmap_gene_count: int,
    upstreams: list[dict],
) -> dict:
    """Assemble the structured health report. Overall status derives from upstreams."""
    if not upstreams:
        overall = "unknown"
    elif all(u.get("reachable") for u in upstreams):
        overall = "degraded" if any(not u.get("ok") for u in upstreams) else "ok"
    elif any(u.get("reachable") for u in upstreams):
        overall = "degraded"
    else:
        overall = "down"
    return {
        "status": overall,
        "version": __version__,
        "n_tools": n_tools,
        "uptime_s": round(uptime_s, 1),
        "depmap_gene_count": depmap_gene_count,
        "upstreams": upstreams,
    }


_STATUS_EMOJI = {"ok": "✅", "degraded": "⚠️", "down": "❌", "unknown": "❔"}


def health_to_markdown(report: dict) -> str:
    """Render a health report as Markdown for the health://status resource."""
    emoji = _STATUS_EMOJI.get(report["status"], "❔")
    lines = [
        "# genesis-bio-mcp — Health",
        f"**Status:** {emoji} {report['status']} | **Version:** {report['version']} | "
        f"**Tools:** {report['n_tools']} | **Uptime:** {report['uptime_s']} s",
        f"**DepMap cache:** {report['depmap_gene_count']} genes loaded",
    ]
    upstreams = report.get("upstreams") or []
    if upstreams:
        reachable = sum(1 for u in upstreams if u.get("reachable"))
        lines += [
            "",
            f"## Upstream reachability ({reachable}/{len(upstreams)} reachable)",
            "| Service | Status | HTTP | Latency |",
            "|---|---|---|---:|",
        ]
        for u in upstreams:
            if not u.get("reachable"):
                state = f"❌ down ({u.get('error', 'error')})"
                http = "—"
            elif u.get("ok"):
                state = "✅ reachable"
                http = str(u.get("status"))
            else:
                state = "⚠️ degraded (5xx)"
                http = str(u.get("status"))
            lines.append(f"| {u['name']} | {state} | {http} | {u.get('latency_ms', '—')} ms |")
    return "\n".join(lines)
