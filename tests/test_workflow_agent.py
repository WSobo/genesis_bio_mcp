"""Unit tests for the workflow agent: tool registry and agent loop."""

from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import AsyncMock, MagicMock, patch

import anthropic
import pytest

from genesis_bio_mcp.workflow_agent import (
    ToolSpec,
    _execute_tool,
    build_tool_registry,
    format_registry_docs,
    run_agent_loop,
)

# ---------------------------------------------------------------------------
# Helpers: build a minimal mock server state
# ---------------------------------------------------------------------------


def _mock_state() -> MagicMock:
    """Return a mock FastMCP state with all expected client attributes."""
    state = MagicMock()
    # Each client method returns an AsyncMock that produces a simple Markdown string
    for attr in (
        "uniprot",
        "open_targets",
        "depmap",
        "gwas",
        "pubchem",
        "chembl",
        "alphafold",
        "string_db",
        "dgidb",
        "clinical_trials",
        "reactome",
        "biogrid",
        "sabdab",
        "gnomad",
        "interpro",
    ):
        client = MagicMock()
        for method in (
            "get_protein",
            "get_association",
            "get_essentiality",
            "get_evidence",
            "get_compounds",
            "get_structure",
            "get_interactome",
            "get_interactions",
            "get_drug_interactions",
            "get_trials",
            "get_pathway_context",
            "get_antibody_structures",
            "get_constraint",
            "get_domains",
        ):
            mock_result = MagicMock()
            mock_result.to_markdown.return_value = f"## Mock {attr}.{method}\n\nData here."
            mock_result.pathways = ["pathway1"]  # for reactome non-empty check
            mock_result.total_partners = 5  # for string_db non-empty check
            async_method = AsyncMock(return_value=mock_result)
            setattr(client, method, async_method)
        setattr(state, attr, client)

    # clinical_trials returns a tuple (trials_list, counts_dict)
    state.clinical_trials.get_trials = AsyncMock(return_value=([], {}))
    # dgidb returns a list of drugs
    state.dgidb.get_drug_interactions = AsyncMock(return_value=[])
    return state


# ---------------------------------------------------------------------------
# Helpers: build fake Anthropic response objects
# ---------------------------------------------------------------------------


def _tool_use_block(name: str, input_kwargs: dict, block_id: str = "tbk_001") -> MagicMock:
    block = MagicMock()
    block.type = "tool_use"
    block.id = block_id
    block.name = name
    block.input = input_kwargs
    return block


def _text_block(text: str) -> MagicMock:
    block = MagicMock()
    block.type = "text"
    block.text = text
    return block


def _make_response(stop_reason: str, content: list) -> MagicMock:
    resp = MagicMock()
    resp.stop_reason = stop_reason
    resp.content = content
    return resp


# ---------------------------------------------------------------------------
# Test 1: Registry construction
# ---------------------------------------------------------------------------


def test_build_tool_registry_has_all_tools():
    """Registry must contain all 15 expected tools and each ToolSpec must be valid."""
    state = _mock_state()
    registry = build_tool_registry(state)

    expected_tools = {
        "resolve_gene",
        "get_protein_info",
        "get_target_disease_association",
        "get_cancer_dependency",
        "get_gwas_evidence",
        "get_compounds",
        "get_chembl_compounds",
        "get_protein_structure",
        "get_protein_interactome",
        "get_biogrid_interactions",
        "get_antibody_structures",
        "get_variant_constraints",
        "get_domain_annotation",
        "get_drug_history",
        "get_pathway_context",
        "get_pathway_members",
        "prioritize_target",
        "compare_targets",
    }
    assert set(registry.keys()) == expected_tools

    # run_biology_workflow is the meta-tool; it must NOT be in the inner registry
    assert "run_biology_workflow" not in registry

    # Every entry must be a ToolSpec with required fields populated
    for name, spec in registry.items():
        assert isinstance(spec, ToolSpec), f"{name} is not a ToolSpec"
        assert spec.name == name, f"ToolSpec.name mismatch for {name}"
        assert spec.description, f"ToolSpec.description empty for {name}"
        assert spec.tool_category, f"ToolSpec.tool_category empty for {name}"
        assert spec.use_when, f"ToolSpec.use_when empty for {name}"
        assert callable(spec.fn), f"ToolSpec.fn not callable for {name}"
        # input_schema must be a valid JSON-schema object
        assert spec.input_schema.get("type") == "object"
        assert "required" in spec.input_schema or "properties" in spec.input_schema


# ---------------------------------------------------------------------------
# Test 2: Single-step agent loop (one tool call then end_turn)
# ---------------------------------------------------------------------------


@pytest.mark.asyncio
async def test_agent_loop_single_tool_call():
    """Agent loop should execute one tool call and return synthesized text."""
    pathway_fn = AsyncMock(return_value="## Pathway Context\n\nBRAF is in MAPK cascade.")

    registry = {
        "get_pathway_context": ToolSpec(
            name="get_pathway_context",
            description="Get pathways",
            input_schema={
                "type": "object",
                "properties": {"gene_symbol": {"type": "string"}},
                "required": ["gene_symbol"],
            },
            tool_category="pathways",
            use_when="Use for pathways",
            fn=pathway_fn,
        )
    }

    # Round 1: Claude requests get_pathway_context
    round1 = _make_response(
        "tool_use",
        [_tool_use_block("get_pathway_context", {"gene_symbol": "BRAF"})],
    )
    # Round 2: Claude returns final answer
    round2 = _make_response(
        "end_turn",
        [_text_block("BRAF participates in the MAPK cascade pathway.")],
    )

    mock_create = AsyncMock(side_effect=[round1, round2])
    mock_messages = MagicMock()
    mock_messages.create = mock_create
    mock_client = MagicMock()
    mock_client.messages = mock_messages

    with patch(
        "genesis_bio_mcp.workflow_agent.anthropic.AsyncAnthropic",
        return_value=mock_client,
    ):
        result = await run_agent_loop("What pathways is BRAF in?", registry)

    assert "MAPK cascade" in result
    pathway_fn.assert_awaited_once_with(gene_symbol="BRAF")
    assert mock_create.await_count == 2


# ---------------------------------------------------------------------------
# Test 3: Multi-step MAPK workflow (three sequential tool calls)
# ---------------------------------------------------------------------------


@pytest.mark.asyncio
async def test_agent_loop_multi_step_mapk_workflow():
    """Agent can chain get_pathway_context → get_drug_history → prioritize_target."""
    pathway_fn = AsyncMock(return_value="## Pathway: MAPK\n\nMAP2K1 is in MAPK cascade.")
    drug_fn = AsyncMock(return_value="## Drug History: MAP2K1\n\nNo approved drugs found.")
    prioritize_fn = AsyncMock(
        return_value="## Target Assessment: MAP2K1\n\nScore: 6.5/10 — Medium."
    )

    registry = {
        "get_pathway_context": ToolSpec(
            name="get_pathway_context",
            description="Get pathways",
            input_schema={
                "type": "object",
                "properties": {"gene_symbol": {"type": "string"}},
                "required": ["gene_symbol"],
            },
            tool_category="pathways",
            use_when="Use for pathways",
            fn=pathway_fn,
        ),
        "get_drug_history": ToolSpec(
            name="get_drug_history",
            description="Get drug history",
            input_schema={
                "type": "object",
                "properties": {"gene_symbol": {"type": "string"}},
                "required": ["gene_symbol"],
            },
            tool_category="druggability",
            use_when="Use for drug landscape",
            fn=drug_fn,
        ),
        "prioritize_target": ToolSpec(
            name="prioritize_target",
            description="Full assessment",
            input_schema={
                "type": "object",
                "properties": {
                    "gene_symbol": {"type": "string"},
                    "indication": {"type": "string"},
                },
                "required": ["gene_symbol", "indication"],
            },
            tool_category="synthesis",
            use_when="Use for full report",
            fn=prioritize_fn,
        ),
    }

    # Step 1: pathway lookup
    step1 = _make_response(
        "tool_use",
        [_tool_use_block("get_pathway_context", {"gene_symbol": "MAP2K1"}, "t1")],
    )
    # Step 2: drug history
    step2 = _make_response(
        "tool_use",
        [_tool_use_block("get_drug_history", {"gene_symbol": "MAP2K1"}, "t2")],
    )
    # Step 3: full prioritization
    step3 = _make_response(
        "tool_use",
        [
            _tool_use_block(
                "prioritize_target",
                {"gene_symbol": "MAP2K1", "indication": "cancer"},
                "t3",
            )
        ],
    )
    # Final synthesis
    final = _make_response(
        "end_turn",
        [_text_block("MAP2K1 is underexplored with no approved drugs. Score: 6.5/10.")],
    )

    mock_create = AsyncMock(side_effect=[step1, step2, step3, final])
    mock_messages = MagicMock()
    mock_messages.create = mock_create
    mock_client = MagicMock()
    mock_client.messages = mock_messages

    with patch(
        "genesis_bio_mcp.workflow_agent.anthropic.AsyncAnthropic",
        return_value=mock_client,
    ):
        result = await run_agent_loop(
            "Find underexplored MAPK targets with no approved drugs", registry
        )

    # All three tools must have been called in order
    assert "MAP2K1" in result
    pathway_fn.assert_awaited_once_with(gene_symbol="MAP2K1")
    drug_fn.assert_awaited_once_with(gene_symbol="MAP2K1")
    prioritize_fn.assert_awaited_once_with(gene_symbol="MAP2K1", indication="cancer")
    assert mock_create.await_count == 4


# ---------------------------------------------------------------------------
# Test 4: Error resilience — tool fn raises, loop continues
# ---------------------------------------------------------------------------


@pytest.mark.asyncio
async def test_agent_loop_tool_error_does_not_crash():
    """If a tool raises, the agent receives an error message and can still respond."""
    broken_fn = AsyncMock(side_effect=RuntimeError("Database connection timeout"))

    registry = {
        "get_protein_info": ToolSpec(
            name="get_protein_info",
            description="Get protein info",
            input_schema={
                "type": "object",
                "properties": {"gene_symbol": {"type": "string"}},
                "required": ["gene_symbol"],
            },
            tool_category="gene_annotation",
            use_when="Use for protein info",
            fn=broken_fn,
        )
    }

    # Claude calls the broken tool
    round1 = _make_response(
        "tool_use",
        [_tool_use_block("get_protein_info", {"gene_symbol": "BRAF"})],
    )
    # Claude handles the error and gives a final answer
    round2 = _make_response(
        "end_turn",
        [_text_block("I was unable to retrieve protein info for BRAF due to a tool error.")],
    )

    mock_create = AsyncMock(side_effect=[round1, round2])
    mock_messages = MagicMock()
    mock_messages.create = mock_create
    mock_client = MagicMock()
    mock_client.messages = mock_messages

    with patch(
        "genesis_bio_mcp.workflow_agent.anthropic.AsyncAnthropic",
        return_value=mock_client,
    ):
        result = await run_agent_loop("Tell me about BRAF protein", registry)

    # Should not raise — must return the final text
    assert isinstance(result, str)
    assert len(result) > 0
    broken_fn.assert_awaited_once()

    # The tool_result passed back to Claude must contain error text, not a crash
    second_call_args = mock_create.call_args_list[1]
    messages_sent = second_call_args.kwargs.get("messages") or second_call_args.args[0]
    # Find the user message with tool_result content
    tool_result_messages = [m for m in messages_sent if isinstance(m.get("content"), list)]
    assert tool_result_messages, "Expected a tool_result message to be sent to Claude"
    tool_result_content = tool_result_messages[-1]["content"][0]["content"]
    assert "error" in tool_result_content.lower() or "Error" in tool_result_content


# ---------------------------------------------------------------------------
# Test 5: _execute_tool with unknown tool name returns error string
# ---------------------------------------------------------------------------


@pytest.mark.asyncio
async def test_execute_tool_unknown_name():
    """_execute_tool must return an error string for unknown tool names, not raise."""
    block = SimpleNamespace(name="nonexistent_tool", input={"gene_symbol": "BRAF"})
    result = await _execute_tool(block, registry={})
    assert "Unknown tool" in result
    assert "nonexistent_tool" in result


# ---------------------------------------------------------------------------
# Test 6: max_iterations safety cap
# ---------------------------------------------------------------------------


@pytest.mark.asyncio
async def test_agent_loop_respects_max_iterations():
    """Agent loop must return a timeout message after max_iterations, not loop forever."""
    # Always return tool_use — never end_turn
    always_tool = _make_response(
        "tool_use",
        [_tool_use_block("get_protein_info", {"gene_symbol": "BRAF"})],
    )
    registry = {
        "get_protein_info": ToolSpec(
            name="get_protein_info",
            description="Get protein info",
            input_schema={
                "type": "object",
                "properties": {"gene_symbol": {"type": "string"}},
                "required": ["gene_symbol"],
            },
            tool_category="gene_annotation",
            use_when="Use for protein info",
            fn=AsyncMock(return_value="## Protein Info"),
        )
    }

    mock_create = AsyncMock(return_value=always_tool)
    mock_messages = MagicMock()
    mock_messages.create = mock_create
    mock_client = MagicMock()
    mock_client.messages = mock_messages

    with patch(
        "genesis_bio_mcp.workflow_agent.anthropic.AsyncAnthropic",
        return_value=mock_client,
    ):
        result = await run_agent_loop("Describe BRAF", registry, max_iterations=3)

    assert "maximum" in result.lower() or "iterations" in result.lower()
    assert mock_create.await_count == 3


# ---------------------------------------------------------------------------
# Test: run_agent_loop returns informative message when ANTHROPIC_API_KEY missing
# ---------------------------------------------------------------------------


@pytest.mark.asyncio
async def test_agent_loop_auth_error_returns_informative_message():
    """AuthenticationError from messages.create must return a human-readable message."""
    registry = {
        "get_protein_info": ToolSpec(
            name="get_protein_info",
            description="Get protein info",
            input_schema={"type": "object", "properties": {}, "required": []},
            tool_category="gene_annotation",
            use_when="Use for protein info",
            fn=AsyncMock(return_value="protein info"),
        ),
    }

    mock_create = AsyncMock(
        side_effect=anthropic.AuthenticationError(
            message="invalid x-api-key",
            response=MagicMock(status_code=401),
            body={"type": "error", "error": {"type": "authentication_error"}},
        )
    )
    mock_messages = MagicMock()
    mock_messages.create = mock_create
    mock_client = MagicMock()
    mock_client.messages = mock_messages

    with patch(
        "genesis_bio_mcp.workflow_agent.anthropic.AsyncAnthropic",
        return_value=mock_client,
    ):
        result = await run_agent_loop("What is BRAF?", registry)

    assert "ANTHROPIC_API_KEY" in result
    assert "claude_desktop_config.json" in result
    # Must not re-raise or propagate
    assert isinstance(result, str)


# ---------------------------------------------------------------------------
# Test 7: format_registry_docs produces valid Markdown with all tools
# ---------------------------------------------------------------------------


def test_format_registry_docs_structure():
    """format_registry_docs must include all tools, grouped by category."""
    state = _mock_state()
    registry = build_tool_registry(state)
    docs = format_registry_docs(registry)

    # Header present
    assert "genesis-bio-mcp Tool Registry" in docs
    assert "18 tools" in docs

    # Every tool name must appear
    for name in registry:
        assert f"`{name}`" in docs, f"Tool '{name}' missing from registry docs"

    # Every category must appear as a section header
    categories = {spec.tool_category for spec in registry.values()}
    for cat in categories:
        assert cat.replace("_", " ").title() in docs, f"Category '{cat}' missing"

    # use_when must appear for each tool
    for spec in registry.values():
        assert spec.use_when[:30] in docs, f"use_when for '{spec.name}' missing"

    # Embedding note must be present
    assert "embedding" in docs.lower()
    assert "use_when" in docs


def test_format_registry_docs_required_inputs():
    """Each tool entry must show its required inputs."""
    state = _mock_state()
    registry = build_tool_registry(state)
    docs = format_registry_docs(registry)

    # prioritize_target requires gene_symbol and indication
    assert "gene_symbol" in docs
    assert "indication" in docs

    # compare_targets requires gene_symbols (array)
    assert "gene_symbols" in docs


# ---------------------------------------------------------------------------
# Test 8: _resolve_symbol helper falls back gracefully on resolution failure
# ---------------------------------------------------------------------------


@pytest.mark.asyncio
async def test_resolve_symbol_fallback_on_exception(monkeypatch):
    """_resolve_symbol must return uppercased input when resolution raises."""
    from genesis_bio_mcp import server

    # Patch _resolve_gene to raise
    async def _broken_resolve(*args, **kwargs):
        raise RuntimeError("UniProt timeout")

    monkeypatch.setattr(server, "_resolve_gene", _broken_resolve)

    # Patch mcp.state.uniprot to anything — it won't be used after the patch
    server.mcp.state = MagicMock()

    symbol, ncbi_id = await server._resolve_symbol("her2")
    assert symbol == "HER2"
    assert ncbi_id is None


@pytest.mark.asyncio
async def test_resolve_symbol_returns_canonical(monkeypatch):
    """_resolve_symbol must return the resolved HGNC symbol and ncbi_gene_id."""
    from genesis_bio_mcp import server
    from genesis_bio_mcp.models import GeneResolution

    async def _mock_resolve(gene_name, *, uniprot_client):
        return GeneResolution(
            hgnc_symbol="ERBB2",
            ncbi_gene_id="2064",
            uniprot_accession="P04626",
            synonyms=["HER2", "NEU"],
            source="uniprot",
        )

    monkeypatch.setattr(server, "_resolve_gene", _mock_resolve)
    server.mcp.state = MagicMock()

    symbol, ncbi_id = await server._resolve_symbol("HER2")
    assert symbol == "ERBB2"
    assert ncbi_id == "2064"
