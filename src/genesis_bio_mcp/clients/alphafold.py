"""AlphaFold Database + RCSB PDB structure client."""

from __future__ import annotations

import logging

import httpx

from genesis_bio_mcp.models import PDBStructure, ProteinStructure

logger = logging.getLogger(__name__)

_ALPHAFOLD_URL = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"

# RCSB search API — find structures by UniProt accession
_RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
_RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"


class AlphaFoldClient:
    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client

    async def get_structure(
        self, gene_symbol: str, uniprot_accession: str | None = None
    ) -> ProteinStructure | None:
        """Return AlphaFold prediction + experimental PDB structures for a gene."""
        accession = uniprot_accession
        if not accession:
            logger.warning("AlphaFold: no UniProt accession provided for %s", gene_symbol)
            return None

        alphafold_plddt, alphafold_url, alphafold_version = await self._fetch_alphafold(accession)
        pdb_structures, total_pdb = await self._fetch_pdb_structures(accession)

        has_ligand = any(s.has_ligand for s in pdb_structures)
        best_res = min(
            (s.resolution_angstrom for s in pdb_structures if s.resolution_angstrom is not None),
            default=None,
        )

        return ProteinStructure(
            gene_symbol=gene_symbol,
            uniprot_accession=accession,
            alphafold_plddt=alphafold_plddt,
            alphafold_model_url=alphafold_url,
            alphafold_version=alphafold_version,
            experimental_structures=pdb_structures,
            total_pdb_structures=total_pdb,
            has_ligand_bound=has_ligand,
            best_resolution=best_res,
        )

    async def _fetch_alphafold(
        self, uniprot_id: str
    ) -> tuple[float | None, str | None, str | None]:
        """Return (mean_plddt, pdb_url, version) from AlphaFold API."""
        try:
            resp = await self._client.get(
                _ALPHAFOLD_URL.format(uniprot_id=uniprot_id),
                timeout=20.0,
            )
            if resp.status_code == 404:
                logger.info("AlphaFold: no model for %s", uniprot_id)
                return None, None, None
            resp.raise_for_status()
            data = resp.json()
            if not data:
                return None, None, None
            entry = data[0]
            plddt = entry.get("meanPlddt") or entry.get("globalMetricValue")
            pdb_url = entry.get("pdbUrl")
            version = entry.get("latestVersion")
            if version:
                version = f"v{version}"
            return (float(plddt) if plddt is not None else None), pdb_url, version
        except Exception as exc:
            logger.warning("AlphaFold API error for %s: %s", uniprot_id, exc)
            return None, None, None

    async def _fetch_pdb_structures(self, uniprot_id: str) -> tuple[list[PDBStructure], int]:
        """Query RCSB for experimental structures linked to a UniProt accession."""
        query = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                    "operator": "exact_match",
                    "value": uniprot_id,
                    "negation": False,
                },
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": 0, "rows": 50},
                "sort": [{"sort_by": "score", "direction": "desc"}],
                "scoring_strategy": "combined",
            },
        }
        try:
            resp = await self._client.post(
                _RCSB_SEARCH_URL,
                json=query,
                timeout=20.0,
            )
            if resp.status_code == 204:
                return [], 0
            resp.raise_for_status()
            result = resp.json()
            hits = result.get("result_set", [])
            total = result.get("total_count", len(hits))
            pdb_ids = [h["identifier"] for h in hits[:20]]
        except Exception as exc:
            logger.warning("RCSB search failed for %s: %s", uniprot_id, exc)
            return [], 0

        structures = []
        for pdb_id in pdb_ids:
            entry = await self._fetch_pdb_entry(pdb_id)
            if entry:
                structures.append(entry)

        structures.sort(key=lambda s: s.resolution_angstrom or 99.0)
        return structures, total

    async def _fetch_pdb_entry(self, pdb_id: str) -> PDBStructure | None:
        """Fetch metadata for a single PDB entry."""
        try:
            resp = await self._client.get(
                _RCSB_ENTRY_URL.format(pdb_id=pdb_id),
                timeout=15.0,
            )
            resp.raise_for_status()
            data = resp.json()

            refine = (data.get("refine") or [{}])[0]
            pdbx = data.get("pdbx_vrpt_summary") or {}
            exptl = (data.get("exptl") or [{}])[0]

            method = exptl.get("method", "UNKNOWN")
            resolution = refine.get("ls_d_res_high") or pdbx.get("PDB_resolution")
            if resolution is not None:
                try:
                    resolution = float(resolution)
                except (ValueError, TypeError):
                    resolution = None

            # Check for bound ligands (non-solvent, non-polymer)
            nonpoly = data.get("rcsb_entry_info", {}).get("nonpolymer_entity_count", 0) or 0
            has_ligand = nonpoly > 0

            # Release year from audit_author or rcsb_accession_info
            year = None
            accession_info = data.get("rcsb_accession_info", {})
            deposit_date = accession_info.get("deposit_date", "")
            if deposit_date and len(deposit_date) >= 4:
                try:
                    year = int(deposit_date[:4])
                except ValueError:
                    pass

            return PDBStructure(
                pdb_id=pdb_id.upper(),
                resolution_angstrom=resolution,
                method=method,
                has_ligand=has_ligand,
                release_year=year,
            )
        except Exception as exc:
            logger.warning("RCSB entry fetch failed for %s: %s", pdb_id, exc)
            return None
