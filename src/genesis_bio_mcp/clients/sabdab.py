"""SAbDab (Structural Antibody Database) client.

SAbDab is a curated database of antibody and nanobody (VHH) crystal structures
from the PDB, maintained by the Oxford Protein Informatics Group (OPIG).

Data access: the summary endpoint returns the full ~20 K row TSV for all structures.
The client caches this locally and filters in-memory to avoid re-downloading on every
tool call.  Cache TTL defaults to 7 days (configurable via GENESIS_SABDAB_CACHE_TTL_SECS).
No API key required.
"""

from __future__ import annotations

import asyncio
import csv
import io
import logging
import time

import httpx

from genesis_bio_mcp.config.settings import settings
from genesis_bio_mcp.models import AntibodyStructure, AntibodyStructures

logger = logging.getLogger(__name__)

_SUMMARY_URL = "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/summary/all/"
_SEMAPHORE = asyncio.Semaphore(2)

# AbNum web API (Chothia scheme) — no install required
_ABNUM_URL = "http://www.bioinf.org.uk/abs/abnum/abnum.cgi"
_ABNUM_SCHEME = "-c"  # Chothia

# RCSB PDB FASTA endpoint for per-entry chain sequences
_RCSB_FASTA_URL = "https://www.rcsb.org/fasta/entry/{pdb}/display"

# CDR position ranges in Chothia heavy-chain numbering (inclusive)
# Insertions (e.g. H52A, H100A) are captured by checking the integer part of the position.
_CHOTHIA_H_CDRS: list[tuple[str, int, int]] = [
    ("vh_cdr1", 26, 32),
    ("vh_cdr2", 52, 58),
    ("vh_cdr3", 95, 102),
]
_CHOTHIA_L_CDRS: list[tuple[str, int, int]] = [
    ("vl_cdr1", 24, 34),
    ("vl_cdr2", 50, 56),
    ("vl_cdr3", 89, 97),
]

# Annotate only the top N structures to keep latency reasonable
_CDR_ANNOTATE_TOP_N = 5

# Fields searched when filtering by query string (case-insensitive substring)
_SEARCH_FIELDS = ("antigen_name", "compound", "antigen_het_name")

# Species strings that indicate a camelid VHH source
_CAMELID_TERMS = ("lama", "camel", "alpaca", "vicugna", "dromedary")


def _parse_tsv(raw: bytes) -> list[dict[str, str]]:
    """Parse SAbDab summary TSV bytes into a list of row dicts."""
    try:
        text = raw.decode("utf-8", errors="replace")
        reader = csv.DictReader(io.StringIO(text), delimiter="\t")
        return [dict(row) for row in reader]
    except Exception as exc:
        logger.warning("SAbDab TSV parse failed: %s", exc)
        return []


def _is_nanobody(row: dict[str, str]) -> bool:
    """Return True if the row represents a VHH/nanobody chain (no light chain)."""
    lchain = row.get("Lchain", "").strip().upper()
    return lchain in ("NA", "NONE", "")


def _parse_resolution(val: str) -> float | None:
    try:
        f = float(val)
        return f if f > 0 else None
    except (ValueError, TypeError):
        return None


def _parse_float(val: str) -> float | None:
    try:
        f = float(val)
        return f if f != 0 else None
    except (ValueError, TypeError):
        return None


def _parse_fasta_chains(fasta_text: str) -> dict[str, str]:
    """Parse RCSB FASTA into {chain_letter: sequence} for all chains in the entry."""
    chains: dict[str, str] = {}
    current_letters: list[str] = []
    seq_lines: list[str] = []

    def _flush() -> None:
        if current_letters and seq_lines:
            seq = "".join(seq_lines)
            for letter in current_letters:
                chains[letter] = seq

    for line in fasta_text.splitlines():
        if line.startswith(">"):
            _flush()
            current_letters = []
            seq_lines = []
            # Header format: >XXXX_N|Chains A, B|description|organism
            parts = line[1:].split("|")
            if len(parts) >= 2 and parts[1].startswith("Chains "):
                current_letters = [c.strip().upper() for c in parts[1][7:].split(",")]
        else:
            seq_lines.append(line.strip())
    _flush()
    return chains


def _parse_abnum(text: str) -> dict[str, str]:
    """Parse AbNum plain-text output into {position_string: residue} dict.

    Each line is ``H26 G`` or ``L31 N`` (chain prefix + position number + optional
    insertion letter, then a space, then the one-letter residue code).
    Returns an empty dict when AbNum returns an error or empty response.
    """
    result: dict[str, str] = {}
    for line in text.strip().splitlines():
        parts = line.split()
        if len(parts) == 2:
            result[parts[0]] = parts[1]
    return result


def _extract_cdrs(
    numbered: dict[str, str],
    chain_prefix: str,
    cdr_ranges: list[tuple[str, int, int]],
) -> dict[str, str | None]:
    """Extract CDR sequences from a numbered antibody chain.

    Args:
        numbered: Output of ``_parse_abnum`` — keys like ``"H26"``, ``"H52A"``.
        chain_prefix: ``"H"`` for heavy, ``"L"`` for light.
        cdr_ranges: List of (field_name, start_pos, end_pos) for Chothia CDRs.

    Returns:
        Dict of field_name → CDR sequence string (or None if no residues found).
    """
    out: dict[str, str | None] = {}
    for field, start, end in cdr_ranges:
        residues: list[tuple[str, str]] = []
        for pos, aa in numbered.items():
            if not pos.startswith(chain_prefix):
                continue
            # Parse numeric part; skip non-numeric suffix for range check
            num_str = "".join(c for c in pos[len(chain_prefix) :] if c.isdigit())
            if not num_str:
                continue
            num = int(num_str)
            if start <= num <= end:
                residues.append((pos, aa))

        if residues:
            # Sort by (numeric_part, insertion_letter) so insertions come in order
            def _sort_key(item: tuple[str, str]) -> tuple[int, str]:
                pos = item[0][len(chain_prefix) :]
                digits = "".join(c for c in pos if c.isdigit())
                letters = "".join(c for c in pos if c.isalpha())
                return (int(digits), letters)

            residues.sort(key=_sort_key)
            out[field] = "".join(aa for _, aa in residues)
        else:
            out[field] = None
    return out


def _row_to_structure(row: dict[str, str]) -> AntibodyStructure:
    nanobody = _is_nanobody(row)
    heavy_species = row.get("heavy_species", "").strip() or None
    light_species = row.get("light_species", "").strip() or None
    affinity_val = row.get("affinity", "None")
    return AntibodyStructure(
        pdb=row.get("pdb", "").strip().upper(),
        is_nanobody=nanobody,
        antigen_name=row.get("antigen_name", "").strip() or None,
        resolution_ang=_parse_resolution(row.get("resolution", "")),
        method=row.get("method", "").strip() or None,
        heavy_species=heavy_species,
        light_species=None if nanobody else light_species,
        heavy_subclass=row.get("heavy_subclass", "").strip() or None,
        light_subclass=None if nanobody else (row.get("light_subclass", "").strip() or None),
        is_engineered=(row.get("engineered", "False").strip().lower() == "true"),
        is_scfv=(row.get("scfv", "False").strip().lower() == "true"),
        affinity_nM=_parse_float(affinity_val) if affinity_val not in ("None", "NA", "") else None,
        compound=row.get("compound", "").strip() or None,
        date_added=row.get("date", "").strip() or None,
        pmid=row.get("pmid", "").strip()
        if row.get("pmid", "").strip() not in ("None", "NA", "")
        else None,
    )


class SAbDabClient:
    """Client for the SAbDab structural antibody database.

    Downloads and caches the full SAbDab summary TSV on first use; subsequent
    queries filter the in-memory database without additional network calls until
    the cache TTL expires.
    """

    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client
        self._db: list[dict[str, str]] | None = None
        self._db_loaded_at: float = 0.0

    async def get_antibody_structures(
        self,
        query: str,
        max_results: int = 20,
    ) -> AntibodyStructures | None:
        """Return antibody/nanobody structures from SAbDab for the given antigen query.

        Searches antigen_name, compound, and antigen_het_name fields (case-insensitive
        substring).  Results sorted by resolution (best first); NA resolution entries
        placed last.

        Args:
            query: Antigen gene symbol or name, e.g. ``"EGFR"`` or ``"epidermal growth factor"``.
            max_results: Maximum number of structures to return (default 20).

        Returns:
            :class:`AntibodyStructures` or ``None`` on download failure.
        """
        await self._ensure_db()
        if self._db is None:
            return None

        q = query.strip().lower()
        hits = [
            row
            for row in self._db
            if any(q in row.get(field, "").lower() for field in _SEARCH_FIELDS)
        ]

        # Dedup by PDB ID — SAbDab has one row per chain pair, so the same
        # crystal can appear twice when both chains match the antigen filter.
        # Keep the best-resolution row per PDB.
        if hits:
            hits.sort(key=lambda r: _parse_resolution(r.get("resolution", "")) or 99.0)
            seen_pdbs: set[str] = set()
            deduped: list[dict[str, str]] = []
            for row in hits:
                pdb = row.get("pdb", "").strip().upper()
                if pdb and pdb not in seen_pdbs:
                    seen_pdbs.add(pdb)
                    deduped.append(row)
            hits = deduped

        if not hits:
            return AntibodyStructures(
                query=query,
                total_structures=0,
                nanobody_count=0,
                fab_count=0,
                structures=[],
            )

        # Sort hits by resolution (best first) keeping rows and structures in sync
        hits_with_structs = sorted(
            ((r, _row_to_structure(r)) for r in hits),
            key=lambda p: (p[1].resolution_ang is None, p[1].resolution_ang or 99.0),
        )

        all_structures = [s for _, s in hits_with_structs]
        top_rows = [r for r, _ in hits_with_structs[:_CDR_ANNOTATE_TOP_N]]
        top_structures = all_structures[:_CDR_ANNOTATE_TOP_N]

        annotated = await self._annotate_cdrs(top_rows, top_structures)
        final_structures = annotated + all_structures[_CDR_ANNOTATE_TOP_N:]

        nanobody_count = sum(1 for s in all_structures if s.is_nanobody)
        fab_count = len(all_structures) - nanobody_count

        return AntibodyStructures(
            query=query,
            total_structures=len(all_structures),
            nanobody_count=nanobody_count,
            fab_count=fab_count,
            structures=final_structures[:max_results],
        )

    # ------------------------------------------------------------------
    # Internal: CDR annotation via RCSB FASTA + AbNum
    # ------------------------------------------------------------------

    async def _annotate_cdrs(
        self,
        rows: list[dict[str, str]],
        structures: list[AntibodyStructure],
    ) -> list[AntibodyStructure]:
        """Annotate structures with CDR sequences by fetching FASTA + running AbNum.

        Fetches RCSB FASTA for each unique PDB concurrently, then submits each
        VH/VL chain sequence to the AbNum web API (Chothia scheme).  Failures are
        silently absorbed — structures are returned unchanged on any error.
        """
        if not structures:
            return structures

        # Fetch all unique FASTA files concurrently
        unique_pdbs = list({r.get("pdb", "").strip().upper() for r in rows if r.get("pdb")})
        fasta_tasks = [self._fetch_fasta(pdb) for pdb in unique_pdbs]
        fasta_results = await asyncio.gather(*fasta_tasks, return_exceptions=True)
        fasta_map: dict[str, dict[str, str]] = {}
        for pdb, result in zip(unique_pdbs, fasta_results):
            if isinstance(result, dict):
                fasta_map[pdb] = result

        annotated: list[AntibodyStructure] = []
        for row, struct in zip(rows, structures):
            pdb = struct.pdb
            chains = fasta_map.get(pdb, {})
            hchain = row.get("Hchain", "").strip().upper()
            lchain = row.get("Lchain", "").strip().upper()

            vh_seq = chains.get(hchain) if hchain else None
            vl_seq = chains.get(lchain) if lchain and not struct.is_nanobody else None

            # Run AbNum for VH and VL concurrently
            vh_task = self._run_abnum(vh_seq, "H") if vh_seq else asyncio.sleep(0)
            vl_task = self._run_abnum(vl_seq, "L") if vl_seq else asyncio.sleep(0)
            vh_cdrs_raw, vl_cdrs_raw = await asyncio.gather(vh_task, vl_task)

            vh_cdrs: dict[str, str | None] = vh_cdrs_raw if isinstance(vh_cdrs_raw, dict) else {}
            vl_cdrs: dict[str, str | None] = vl_cdrs_raw if isinstance(vl_cdrs_raw, dict) else {}

            annotated.append(
                struct.model_copy(
                    update={
                        "vh_cdr1": vh_cdrs.get("vh_cdr1"),
                        "vh_cdr2": vh_cdrs.get("vh_cdr2"),
                        "vh_cdr3": vh_cdrs.get("vh_cdr3"),
                        "vl_cdr1": vl_cdrs.get("vl_cdr1"),
                        "vl_cdr2": vl_cdrs.get("vl_cdr2"),
                        "vl_cdr3": vl_cdrs.get("vl_cdr3"),
                    }
                )
            )
        return annotated

    async def _fetch_fasta(self, pdb: str) -> dict[str, str]:
        """Fetch RCSB FASTA for *pdb* and return {chain_letter: sequence}."""
        url = _RCSB_FASTA_URL.format(pdb=pdb)
        try:
            resp = await self._client.get(url, timeout=20.0)
            resp.raise_for_status()
            return _parse_fasta_chains(resp.text)
        except Exception as exc:
            logger.warning("SAbDab: RCSB FASTA fetch failed for %s: %s", pdb, exc)
            return {}

    async def _run_abnum(self, sequence: str, chain_prefix: str) -> dict[str, str | None]:
        """Submit *sequence* to AbNum (Chothia) and return extracted CDR dict."""
        cdr_ranges = _CHOTHIA_H_CDRS if chain_prefix == "H" else _CHOTHIA_L_CDRS
        try:
            resp = await self._client.get(
                _ABNUM_URL,
                params={"plain": "1", "scheme": _ABNUM_SCHEME, "aaseq": sequence},
                timeout=20.0,
            )
            resp.raise_for_status()
            numbered = _parse_abnum(resp.text)
            if not numbered:
                return {field: None for field, _, _ in cdr_ranges}
            return _extract_cdrs(numbered, chain_prefix, cdr_ranges)
        except Exception as exc:
            logger.warning("SAbDab: AbNum failed for chain %s: %s", chain_prefix, exc)
            return {field: None for field, _, _ in cdr_ranges}

    # ------------------------------------------------------------------
    # Internal: cache management
    # ------------------------------------------------------------------

    async def _ensure_db(self) -> None:
        """Load (or refresh) the full SAbDab TSV into memory."""
        now = time.time()
        ttl = settings.sabdab_cache_ttl_secs

        # In-memory still fresh
        if self._db is not None and (now - self._db_loaded_at) < ttl:
            logger.debug("SAbDab: using in-memory cache (%d rows)", len(self._db))
            return

        # Try disk cache
        cache_path = settings.sabdab_cache_path
        if cache_path.exists():
            age = now - cache_path.stat().st_mtime
            if age < ttl:
                logger.debug("SAbDab: loading from disk cache at %s", cache_path)
                try:
                    raw = cache_path.read_bytes()
                    self._db = _parse_tsv(raw)
                    self._db_loaded_at = now
                    return
                except Exception as exc:
                    logger.warning("SAbDab: disk cache read failed: %s", exc)

        # Download fresh
        async with _SEMAPHORE:
            await self._download_and_cache(cache_path)

    async def _download_and_cache(self, cache_path) -> None:
        """Download the SAbDab summary TSV and persist to disk."""
        logger.info("SAbDab: downloading fresh summary TSV from %s", _SUMMARY_URL)
        try:
            resp = await self._client.get(_SUMMARY_URL, timeout=60.0)
            resp.raise_for_status()
            raw = resp.content
        except Exception as exc:
            logger.warning("SAbDab: download failed: %s", exc)
            # Keep stale cache if available
            if cache_path.exists():
                logger.warning("SAbDab: falling back to stale disk cache")
                try:
                    self._db = _parse_tsv(cache_path.read_bytes())
                    self._db_loaded_at = time.time()
                except Exception:
                    pass
            return

        # Persist to disk
        try:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            cache_path.write_bytes(raw)
        except Exception as exc:
            logger.warning("SAbDab: failed to write disk cache: %s", exc)

        self._db = _parse_tsv(raw)
        self._db_loaded_at = time.time()
        logger.info("SAbDab: loaded %d rows", len(self._db or []))
