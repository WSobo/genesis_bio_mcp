"""InterPro domain annotation client.

InterPro is a database of protein families, domains, and functional sites that
integrates predictive models from multiple member databases (Pfam, SMART, CDD,
PRINTS, PANTHER, etc.).

This client fetches InterPro domain annotations for a protein given its UniProt
accession.  Knowing domain boundaries is critical for protein engineering: avoid
engineering within conserved domain cores; focus mutagenesis on loop regions or
linkers between domains.

No API key required.  Uses the EBI InterPro REST API v4.
"""

from __future__ import annotations

import asyncio
import logging

import httpx

from genesis_bio_mcp.models import DomainAnnotation, DomainAnnotations

logger = logging.getLogger(__name__)

_INTERPRO_BASE = "https://www.ebi.ac.uk/interpro/api"
_SEMAPHORE = asyncio.Semaphore(3)


class InterProClient:
    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client
        self._cache: dict[str, DomainAnnotations] = {}

    async def get_domains(
        self,
        gene_symbol: str,
        uniprot_accession: str,
        page_size: int = 30,
    ) -> DomainAnnotations | None:
        """Return InterPro domain annotations for *uniprot_accession*.

        Fetches only InterPro-integrated entries (not member-database-only entries)
        to avoid redundancy.  Results include domain name, accession, position,
        member database cross-references, and GO terms.

        Args:
            gene_symbol: HGNC symbol (used for display and caching key).
            uniprot_accession: Primary UniProt accession, e.g. ``'P15056'``.
            page_size: Maximum number of domain entries to retrieve (default 30).

        Returns:
            :class:`DomainAnnotations` or ``None`` on error.
        """
        key = uniprot_accession.upper()
        if key in self._cache:
            logger.debug("InterPro cache hit: %s", key)
            return self._cache[key]

        async with _SEMAPHORE:
            result = await self._fetch(gene_symbol, key, page_size)

        if result is not None:
            self._cache[key] = result
        return result

    async def _fetch(
        self, gene_symbol: str, accession: str, page_size: int
    ) -> DomainAnnotations | None:
        url = f"{_INTERPRO_BASE}/entry/InterPro/protein/UniProt/{accession}/"
        try:
            resp = await self._client.get(
                url,
                params={"page_size": page_size},
                timeout=20.0,
            )
            if resp.status_code == 404:
                logger.info("InterPro: no entries for accession '%s'", accession)
                return DomainAnnotations(
                    gene_symbol=gene_symbol,
                    uniprot_accession=accession,
                    total_entries=0,
                    domains=[],
                )
            resp.raise_for_status()
            data = resp.json()
        except Exception as exc:
            logger.warning("InterPro fetch failed for %s (%s): %s", gene_symbol, accession, exc)
            return None

        entries = []
        for result in data.get("results", []):
            meta = result.get("metadata", {})
            proteins = result.get("proteins", [])
            locations = proteins[0].get("entry_protein_locations", []) if proteins else []

            # Collect residue positions from all location fragments
            positions: list[tuple[int, int]] = []
            for loc in locations:
                for frag in loc.get("fragments", []):
                    start = frag.get("start")
                    end = frag.get("end")
                    if start is not None and end is not None:
                        positions.append((int(start), int(end)))

            # Member database cross-refs
            member_dbs: dict[str, list[str]] = {}
            for db_name, db_entries in (meta.get("member_databases") or {}).items():
                member_dbs[db_name] = list(db_entries.keys())

            # GO terms
            go_terms = [f"{g['identifier']} {g['name']}" for g in (meta.get("go_terms") or [])]

            entries.append(
                DomainAnnotation(
                    interpro_accession=meta.get("accession", ""),
                    name=meta.get("name", ""),
                    entry_type=meta.get("type", ""),
                    positions=positions,
                    member_databases=member_dbs,
                    go_terms=go_terms,
                )
            )

        # Sort by start position of first fragment
        entries.sort(key=lambda e: e.positions[0][0] if e.positions else 9999)

        return DomainAnnotations(
            gene_symbol=gene_symbol,
            uniprot_accession=accession,
            total_entries=data.get("count", len(entries)),
            domains=entries,
        )
