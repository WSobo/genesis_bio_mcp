"""Cross-database gene symbol resolution.

Resolution priority chain:
1. UniProt gene_exact search → HGNC symbol + UniProt accession
2. NCBI E-utils esearch → NCBI Gene ID
3. UniProt gene_synonym search for alias resolution
"""

from __future__ import annotations

import logging
import os
from typing import TYPE_CHECKING

import httpx

from genesis_bio_mcp.models import GeneResolution

if TYPE_CHECKING:
    from genesis_bio_mcp.clients.uniprot import UniProtClient

logger = logging.getLogger(__name__)

_NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
_NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "genesis-bio-mcp@example.com")


async def resolve_gene(
    gene_name: str,
    *,
    uniprot_client: UniProtClient,
    http_client: httpx.AsyncClient | None = None,
) -> GeneResolution:
    """Resolve a gene name or alias to canonical identifiers.

    Falls back gracefully — always returns a GeneResolution even if
    some fields cannot be populated.
    """
    symbol = gene_name.strip().upper()

    # Step 1: UniProt primary search
    uniprot_entry = await uniprot_client._search(symbol, reviewed_only=True)
    if uniprot_entry is None:
        uniprot_entry = await uniprot_client._search(symbol, reviewed_only=False)

    hgnc_symbol = symbol
    uniprot_accession: str | None = None
    synonyms: list[str] = []
    source = "input"

    if uniprot_entry:
        hgnc_symbol, uniprot_accession, synonyms = _extract_gene_info(uniprot_entry)
        source = "uniprot"
    else:
        # Step 3: Try synonym search
        synonym_entry = await uniprot_client.search_by_synonym(symbol)
        if synonym_entry:
            hgnc_symbol, uniprot_accession, synonyms = _extract_gene_info(synonym_entry)
            source = "uniprot_synonym"

    # Step 2: NCBI E-utils for gene ID
    ncbi_gene_id: str | None = None
    client = http_client or uniprot_client._client
    ncbi_gene_id = await _fetch_ncbi_gene_id(hgnc_symbol, client)

    return GeneResolution(
        hgnc_symbol=hgnc_symbol,
        hgnc_id=None,  # Would require HGNC REST API; omitted for simplicity
        ncbi_gene_id=ncbi_gene_id,
        uniprot_accession=uniprot_accession,
        synonyms=synonyms,
        source=source,
    )


def _extract_gene_info(entry: dict) -> tuple[str, str | None, list[str]]:
    accession = entry.get("primaryAccession")
    genes = entry.get("genes", [])
    symbol = ""
    synonyms: list[str] = []

    if genes:
        gn = genes[0]
        symbol = gn.get("geneName", {}).get("value", "")
        for s in gn.get("synonyms", []):
            v = s.get("value")
            if v and v != symbol:
                synonyms.append(v)

    return symbol.upper() if symbol else "", accession, synonyms


async def _fetch_ncbi_gene_id(symbol: str, client: httpx.AsyncClient) -> str | None:
    params = {
        "db": "gene",
        "term": f"{symbol}[Gene Name] AND Homo sapiens[Organism]",
        "retmode": "json",
        "retmax": "1",
        "email": _NCBI_EMAIL,
    }
    try:
        resp = await client.get(_NCBI_ESEARCH, params=params, timeout=10.0)
        resp.raise_for_status()
        id_list = resp.json().get("esearchresult", {}).get("idlist", [])
        return id_list[0] if id_list else None
    except Exception as exc:
        logger.debug("NCBI E-utils gene ID lookup failed for '%s': %s", symbol, exc)
        return None
