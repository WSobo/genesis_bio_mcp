"""Tests for the health://status resource helpers."""

import httpx
import respx

from genesis_bio_mcp.tools.health import build_health_report, check_upstreams, health_to_markdown


def _ok(name: str, status: int = 200) -> dict:
    return {"name": name, "reachable": True, "status": status, "latency_ms": 50, "ok": status < 500}


def _down(name: str) -> dict:
    return {
        "name": name,
        "reachable": False,
        "status": None,
        "latency_ms": 8000,
        "ok": False,
        "error": "ConnectError",
    }


def test_build_health_report_status_logic():
    base = dict(n_tools=36, uptime_s=12.3, depmap_gene_count=18000)
    assert build_health_report(**base, upstreams=[_ok("A"), _ok("B")])["status"] == "ok"
    # reachable but 5xx → degraded
    assert build_health_report(**base, upstreams=[_ok("A"), _ok("B", 503)])["status"] == "degraded"
    # some reachable, some down → degraded
    assert build_health_report(**base, upstreams=[_ok("A"), _down("B")])["status"] == "degraded"
    # none reachable → down
    assert build_health_report(**base, upstreams=[_down("A"), _down("B")])["status"] == "down"
    # no probes → unknown
    assert build_health_report(**base, upstreams=[])["status"] == "unknown"

    report = build_health_report(**base, upstreams=[_ok("A")])
    assert report["n_tools"] == 36
    assert report["depmap_gene_count"] == 18000
    assert report["version"]  # non-empty version string


def test_health_to_markdown_renders():
    report = build_health_report(
        n_tools=36,
        uptime_s=12.3,
        depmap_gene_count=18000,
        upstreams=[_ok("UniProt"), _ok("OT", 400), _down("DownSvc")],
    )
    md = health_to_markdown(report)
    assert "genesis-bio-mcp — Health" in md
    assert "**Tools:** 36" in md
    assert "18000 genes loaded" in md
    assert "UniProt" in md and "✅ reachable" in md
    assert "OT" in md  # HTTP 400 is still reachable (< 500)
    assert report["status"] == "degraded"  # one upstream down → overall degraded
    assert "DownSvc" in md and "down" in md


@respx.mock
async def test_check_upstreams_reports_reachability(http_client):
    # Specific (down) route first so it takes precedence over the catch-all.
    respx.get(url__regex=r"string-db\.org").mock(side_effect=httpx.ConnectError("boom"))
    respx.get(url__regex=r".+").mock(return_value=httpx.Response(200, json={}))

    ups = await check_upstreams(http_client)
    by_name = {u["name"]: u for u in ups}
    assert by_name["UniProt"]["reachable"] is True
    assert by_name["UniProt"]["ok"] is True
    assert by_name["STRING"]["reachable"] is False
    assert by_name["STRING"]["error"] == "ConnectError"
