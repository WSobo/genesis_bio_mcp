"""Fetch the human kinome (targets + sequences) from UniProt for corpus ingestion.

The corpus's target family is defined by a UniProt query — reviewed human proteins carrying
the Kinase keyword (KW-0418). One paged search returns accession + gene + sequence +
protein name per target; no per-protein round-trips. Public REST, no API key.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import httpx

logger = logging.getLogger(__name__)

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
# Reviewed (Swiss-Prot) human proteins with the Kinase keyword. This IS the corpus's
# definition of "the kinome" — recorded so the build is reproducible.
DEFAULT_QUERY = "(keyword:KW-0418) AND (organism_id:9606) AND (reviewed:true)"
_FIELDS = "accession,gene_primary,protein_name,sequence,organism_name"


@dataclass
class TargetRecord:
    """A target protein to index: UniProt-keyed, with the sequence to embed."""

    uniprot_accession: str
    gene_symbol: str | None
    protein_name: str | None
    organism: str | None
    sequence: str


def _parse_entry(entry: dict) -> TargetRecord | None:
    """Map one UniProtKB JSON entry to a TargetRecord, or None if it lacks accession/sequence."""
    acc = entry.get("primaryAccession")
    seq = (entry.get("sequence") or {}).get("value")
    if not acc or not seq:
        return None
    genes = entry.get("genes") or []
    gene = (genes[0].get("geneName") or {}).get("value") if genes else None
    pname = (
        (((entry.get("proteinDescription") or {}).get("recommendedName") or {}).get("fullName"))
        or {}
    ).get("value")
    organism = (entry.get("organism") or {}).get("scientificName")
    return TargetRecord(
        uniprot_accession=acc,
        gene_symbol=gene,
        protein_name=pname,
        organism=organism,
        sequence=seq,
    )


def _next_link(link_header: str | None) -> str | None:
    """Extract the rel="next" URL from a UniProt Link header for cursor pagination."""
    if not link_header:
        return None
    for part in link_header.split(","):
        section = part.split(";")
        if len(section) < 2:
            continue
        url = section[0].strip().lstrip("<").rstrip(">")
        if 'rel="next"' in section[1]:
            return url
    return None


async def fetch_kinome_targets(
    client: httpx.AsyncClient,
    *,
    query: str = DEFAULT_QUERY,
    page_size: int = 500,
    max_targets: int | None = None,
) -> list[TargetRecord]:
    """Fetch kinome target records from UniProt, following cursor pagination.

    ``max_targets`` caps the result (handy for a quick partial build); ``None`` fetches all.
    """
    records: list[TargetRecord] = []
    next_url: str | None = None
    params: dict | None = {
        "query": query,
        "fields": _FIELDS,
        "format": "json",
        "size": str(page_size),
    }
    while True:
        resp = await client.get(next_url or UNIPROT_SEARCH, params=None if next_url else params)
        resp.raise_for_status()
        for entry in resp.json().get("results", []):
            rec = _parse_entry(entry)
            if rec is not None:
                records.append(rec)
                if max_targets is not None and len(records) >= max_targets:
                    logger.info("Fetched %d targets (capped).", len(records))
                    return records
        next_url = _next_link(resp.headers.get("link"))
        if not next_url:
            break
    logger.info("Fetched %d kinome targets from UniProt.", len(records))
    return records
