"""AlphaFold Database + RCSB PDB structure client."""

from __future__ import annotations

import asyncio
import logging

import httpx

from genesis_bio_mcp.models import (
    ConfidenceRegion,
    PDBStructure,
    ProteinStructure,
    StructureConfidence,
)

logger = logging.getLogger(__name__)

# Cap concurrent EBI AlphaFold + RCSB fetches.
_SEMAPHORE = asyncio.Semaphore(3)

_ALPHAFOLD_URL = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"

# RCSB search API — find structures by UniProt accession
_RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
_RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"


class AlphaFoldClient:
    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client
        # Session-scoped cache keyed by (gene_symbol, uniprot_accession).
        # Structure data is stable within a session; PDB deposits don't change mid-run.
        self._cache: dict[tuple[str, str | None], ProteinStructure | None] = {}
        # Separate cache for the heavier per-residue confidence profile (downloads
        # and parses the full model PDB), keyed by (gene_symbol, accession).
        self._confidence_cache: dict[tuple[str, str], StructureConfidence | None] = {}

    async def get_structure(
        self, gene_symbol: str, uniprot_accession: str | None = None
    ) -> ProteinStructure | None:
        """Return AlphaFold prediction + experimental PDB structures for a gene."""
        accession = uniprot_accession
        if not accession:
            logger.warning("AlphaFold: no UniProt accession provided for %s", gene_symbol)
            return None

        cache_key = (gene_symbol.upper(), accession)
        if cache_key in self._cache:
            logger.debug("AlphaFold cache hit: %s", gene_symbol)
            return self._cache[cache_key]

        alphafold_plddt, alphafold_url, alphafold_version = await self._fetch_alphafold(accession)
        pdb_structures, total_pdb = await self._fetch_pdb_structures(accession)

        has_ligand = any(s.has_ligand for s in pdb_structures)
        best_res = min(
            (s.resolution_angstrom for s in pdb_structures if s.resolution_angstrom is not None),
            default=None,
        )

        result = ProteinStructure(
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
        self._cache[cache_key] = result
        return result

    async def get_confidence(
        self, gene_symbol: str, uniprot_accession: str | None = None
    ) -> StructureConfidence | None:
        """Return a per-residue pLDDT confidence profile from the AlphaFold model.

        Downloads the AlphaFold model PDB and parses per-residue pLDDT from the
        CA-atom B-factor column, then summarizes confidence bands and contiguous
        low-confidence (likely flexible/disordered) regions.
        """
        if not uniprot_accession:
            logger.warning("AlphaFold confidence: no UniProt accession for %s", gene_symbol)
            return None

        cache_key = (gene_symbol.upper(), uniprot_accession)
        if cache_key in self._confidence_cache:
            logger.debug("AlphaFold confidence cache hit: %s", gene_symbol)
            return self._confidence_cache[cache_key]

        _, pdb_url, version = await self._fetch_alphafold(uniprot_accession)
        if not pdb_url:
            logger.info("AlphaFold confidence: no model for %s", uniprot_accession)
            self._confidence_cache[cache_key] = None
            return None

        try:
            async with _SEMAPHORE:
                resp = await self._client.get(pdb_url, timeout=25.0)
            resp.raise_for_status()
            plddts = _parse_plddt_from_pdb(resp.text)
        except Exception as exc:
            logger.warning("AlphaFold model download failed for %s: %s", uniprot_accession, exc)
            return None

        if not plddts:
            logger.warning("AlphaFold confidence: no CA pLDDT parsed for %s", uniprot_accession)
            self._confidence_cache[cache_key] = None
            return None

        result = _build_confidence(gene_symbol, uniprot_accession, version, plddts)
        self._confidence_cache[cache_key] = result
        return result

    async def get_model_pdb(self, uniprot_accession: str) -> str | None:
        """Return the raw AlphaFold model PDB text for a UniProt accession.

        Used by structural-search tooling (e.g. Foldseek) that needs the actual
        coordinate file rather than a parsed summary. Returns ``None`` when the
        protein has no AlphaFold model or the download fails.
        """
        if not uniprot_accession:
            return None
        _, pdb_url, _ = await self._fetch_alphafold(uniprot_accession)
        if not pdb_url:
            logger.info("AlphaFold: no model PDB URL for %s", uniprot_accession)
            return None
        try:
            async with _SEMAPHORE:
                resp = await self._client.get(pdb_url, timeout=25.0)
            resp.raise_for_status()
            return resp.text
        except Exception as exc:
            logger.warning("AlphaFold model download failed for %s: %s", uniprot_accession, exc)
            return None

    async def _fetch_alphafold(
        self, uniprot_id: str
    ) -> tuple[float | None, str | None, str | None]:
        """Return (mean_plddt, pdb_url, version) from AlphaFold API."""
        try:
            async with _SEMAPHORE:
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
            async with _SEMAPHORE:
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
            async with _SEMAPHORE:
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


def _parse_plddt_from_pdb(pdb_text: str) -> list[tuple[int, float]]:
    """Extract ordered (residue_number, pLDDT) from CA-atom B-factors of an AlphaFold PDB.

    AlphaFold stores the per-residue pLDDT in the temperature-factor (B-factor)
    column. Reads one value per residue (the CA atom), preserving sequence order.
    """
    out: list[tuple[int, float]] = []
    seen: set[int] = set()
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != "CA":
            continue
        try:
            res_num = int(line[22:26])
            plddt = float(line[60:66])
        except (ValueError, IndexError):
            continue
        if res_num in seen:
            continue
        seen.add(res_num)
        out.append((res_num, plddt))
    return out


def _find_low_confidence_regions(
    plddts: list[tuple[int, float]], threshold: float, min_length: int
) -> list[ConfidenceRegion]:
    """Group consecutive residues with pLDDT < threshold into regions (≥ min_length)."""
    regions: list[ConfidenceRegion] = []
    run: list[tuple[int, float]] = []

    def flush() -> None:
        nonlocal run
        if len(run) >= min_length:
            vals = [v for _, v in run]
            regions.append(
                ConfidenceRegion(
                    start=run[0][0],
                    end=run[-1][0],
                    length=len(run),
                    mean_plddt=round(sum(vals) / len(vals), 1),
                )
            )
        run = []

    for res_num, plddt in plddts:
        if plddt < threshold:
            if run and res_num != run[-1][0] + 1:
                flush()
            run.append((res_num, plddt))
        else:
            flush()
    flush()
    return regions


def _build_confidence(
    gene_symbol: str,
    uniprot_accession: str,
    version: str | None,
    plddts: list[tuple[int, float]],
) -> StructureConfidence:
    """Summarize per-residue pLDDT into band fractions and low-confidence regions."""
    values = [p for _, p in plddts]
    n = len(values)
    very_high = sum(1 for v in values if v >= 90)
    confident = sum(1 for v in values if 70 <= v < 90)
    low = sum(1 for v in values if 50 <= v < 70)
    very_low = sum(1 for v in values if v < 50)

    return StructureConfidence(
        gene_symbol=gene_symbol,
        uniprot_accession=uniprot_accession,
        model_version=version,
        residue_count=n,
        mean_plddt=round(sum(values) / n, 2),
        pct_very_high=round(100 * very_high / n, 1),
        pct_confident=round(100 * confident / n, 1),
        pct_low=round(100 * low / n, 1),
        pct_very_low=round(100 * very_low / n, 1),
        low_confidence_regions=_find_low_confidence_regions(plddts, threshold=70.0, min_length=3),
    )
