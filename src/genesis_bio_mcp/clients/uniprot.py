"""UniProt REST API client for protein annotation retrieval."""

from __future__ import annotations

import logging

import httpx

from genesis_bio_mcp.models import KnownVariant, ProteinInfo

logger = logging.getLogger(__name__)

_BASE_URL = "https://rest.uniprot.org/uniprotkb"
_FIELDS = (
    "accession,gene_names,protein_name,organism_name,"
    "cc_function,cc_subcellular_location,cc_disease,"
    "xref_reactome,xref_pdb,ft_variant,ft_disulfid"
)


class UniProtClient:
    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client
        # Session-scoped cache: UniProt Swiss-Prot entries are stable within a run.
        self._cache: dict[str, ProteinInfo | None] = {}
        # Separate FASTA cache keyed by (accession, start, end) — sequences are
        # large enough that re-fetching after a slice would be wasteful.
        self._fasta_cache: dict[tuple[str, int | None, int | None], tuple[str, str, str]] = {}

    async def get_protein(self, gene_symbol: str) -> ProteinInfo | None:
        """Return Swiss-Prot annotation for a human gene symbol, or None if not found."""
        symbol = gene_symbol.strip().upper()
        if symbol in self._cache:
            logger.debug("UniProt cache hit: %s", symbol)
            return self._cache[symbol]
        data = await self._search(symbol, reviewed_only=True)
        if data is None:
            data = await self._search(symbol, reviewed_only=False)
        result = _parse_entry(data, symbol) if data is not None else None
        self._cache[symbol] = result
        return result

    async def get_sequence(
        self,
        accession: str,
        start: int | None = None,
        end: int | None = None,
    ) -> tuple[str, str, str] | None:
        """Fetch the FASTA sequence for a UniProt accession.

        Args:
            accession: Primary UniProt accession, e.g. ``"P15056"``.
            start: Optional 1-indexed inclusive slice start.
            end: Optional 1-indexed inclusive slice end (must be ≥ start).

        Returns:
            Tuple of ``(sequence, organism, description)`` on success, or
            ``None`` if the FASTA is missing or unreachable. ``organism`` and
            ``description`` are parsed from the FASTA header; both are best-
            effort ("" if parsing fails).
        """
        acc = accession.strip()
        if not acc:
            return None
        cache_key = (acc, start, end)
        if cache_key in self._fasta_cache:
            logger.debug("UniProt FASTA cache hit: %s", cache_key)
            return self._fasta_cache[cache_key]
        url = f"{_BASE_URL}/{acc}.fasta"
        try:
            resp = await self._client.get(url, timeout=20.0)
            resp.raise_for_status()
        except Exception as exc:
            logger.warning("UniProt FASTA fetch failed for %s: %s", acc, exc)
            return None
        parsed = _parse_fasta(resp.text)
        if parsed is None:
            return None
        full_seq, organism, description = parsed
        if start is not None or end is not None:
            s = start - 1 if start is not None else 0
            e = end if end is not None else len(full_seq)
            if s < 0 or e > len(full_seq) or s >= e:
                logger.warning(
                    "UniProt FASTA slice out of range: accession=%s start=%s end=%s length=%d",
                    acc,
                    start,
                    end,
                    len(full_seq),
                )
                return None
            seq = full_seq[s:e]
        else:
            seq = full_seq
        result = (seq, organism, description)
        self._fasta_cache[cache_key] = result
        return result

    async def search_by_synonym(self, synonym: str) -> dict | None:
        """Search UniProt by gene synonym/alias and return the raw first result."""
        symbol = synonym.strip().upper()
        params = {
            "query": f"gene_synonym:{symbol} AND organism_id:9606 AND reviewed:true",
            "fields": _FIELDS,
            "format": "json",
            "size": "1",
        }
        try:
            resp = await self._client.get(f"{_BASE_URL}/search", params=params)
            resp.raise_for_status()
            results = resp.json().get("results", [])
            return results[0] if results else None
        except Exception as exc:
            logger.warning("UniProt synonym search failed for '%s': %s", synonym, exc)
            return None

    async def _search(self, symbol: str, *, reviewed_only: bool) -> dict | None:
        query = f"gene_exact:{symbol} AND organism_id:9606"
        if reviewed_only:
            query += " AND reviewed:true"
        params = {"query": query, "fields": _FIELDS, "format": "json", "size": "1"}
        try:
            resp = await self._client.get(f"{_BASE_URL}/search", params=params)
            resp.raise_for_status()
            results = resp.json().get("results", [])
            return results[0] if results else None
        except Exception as exc:
            logger.warning(
                "UniProt search failed for '%s' (reviewed=%s): %s",
                symbol,
                reviewed_only,
                exc,
            )
            return None


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------


def _parse_entry(entry: dict, gene_symbol: str) -> ProteinInfo:
    accession = entry.get("primaryAccession", "")
    reviewed = entry.get("entryType", "") == "UniProtKB reviewed (Swiss-Prot)"

    # Protein name
    pn = entry.get("proteinDescription", {})
    rec = pn.get("recommendedName", {})
    protein_name = rec.get("fullName", {}).get("value", "") if rec else ""
    if not protein_name:
        # Fall back to submitted name
        submitted = pn.get("submissionNames", [])
        if submitted:
            protein_name = submitted[0].get("fullName", {}).get("value", "")

    organism = entry.get("organism", {}).get("scientificName", "Homo sapiens")

    # Gene symbol from entry
    genes = entry.get("genes", [])
    canonical_symbol = gene_symbol
    if genes:
        gn = genes[0].get("geneName", {})
        canonical_symbol = gn.get("value", gene_symbol)

    # Comments
    comments = entry.get("comments", [])
    function_summary = ""
    subcellular_locations: list[str] = []
    disease_associations: list[str] = []

    for comment in comments:
        ctype = comment.get("commentType", "")
        if ctype == "FUNCTION":
            texts = comment.get("texts", [])
            if texts:
                function_summary = texts[0].get("value", "")
        elif ctype == "SUBCELLULAR LOCATION":
            for sl in comment.get("subcellularLocations", []):
                loc = sl.get("location", {}).get("value")
                if loc:
                    subcellular_locations.append(loc)
        elif ctype == "DISEASE":
            disease = comment.get("disease", {})
            dname = disease.get("diseaseName")
            if dname:
                disease_associations.append(dname)

    # Cross-references
    xrefs = entry.get("uniProtKBCrossReferences", [])
    pdb_structures: list[str] = []
    pathways: list[str] = []
    for xref in xrefs:
        db = xref.get("database", "")
        xid = xref.get("id", "")
        if db == "PDB":
            pdb_structures.append(xid)
        elif db == "Reactome":
            props = xref.get("properties", [])
            pathway_name = next(
                (p.get("value") for p in props if p.get("key") == "PathwayName"), xid
            )
            pathways.append(pathway_name)

    # Variants + disulfide-bond features
    features = entry.get("features", [])
    known_variants: list[KnownVariant] = []
    disulfide_bond_positions: list[int] = []
    for feat in features:
        ftype = feat.get("type")
        if ftype == "Natural variant":
            loc = feat.get("location", {})
            pos = str(loc.get("start", {}).get("value", ""))
            orig = feat.get("alternativeSequence", {}).get("originalSequence", "")
            alts = feat.get("alternativeSequence", {}).get("alternativeSequences", [])
            alt = alts[0] if alts else ""
            desc = feat.get("description", "")
            disease = desc if "In" in desc else None
            known_variants.append(
                KnownVariant(position=pos, original=orig, variant=alt, disease=disease)
            )
        elif ftype == "Disulfide bond":
            # UniProt reports disulfide bonds as features with both start and
            # end positions pointing at the two bonded Cys residues.
            loc = feat.get("location", {})
            start_val = loc.get("start", {}).get("value")
            end_val = loc.get("end", {}).get("value")
            for v in (start_val, end_val):
                if isinstance(v, int):
                    disulfide_bond_positions.append(v)

    return ProteinInfo(
        uniprot_accession=accession,
        gene_symbol=canonical_symbol,
        protein_name=protein_name,
        organism=organism,
        function_summary=function_summary,
        subcellular_locations=list(dict.fromkeys(subcellular_locations)),
        pathways=list(dict.fromkeys(pathways)),
        disease_associations=disease_associations,
        pdb_structures=pdb_structures[:20],
        known_variants=known_variants[:10],
        disulfide_bond_positions=sorted(set(disulfide_bond_positions)),
        reviewed=reviewed,
    )


def _parse_fasta(text: str) -> tuple[str, str, str] | None:
    """Parse a single-entry UniProt FASTA into (sequence, organism, description).

    UniProt headers follow the format::

        >sp|P15056|BRAF_HUMAN Serine/threonine-protein kinase B-raf OS=Homo sapiens OX=9606 GN=BRAF PE=1 SV=4

    We extract the free-text description (between the third ``|`` and the
    first ``OS=``) and the organism (between ``OS=`` and ``OX=``). Both are
    best-effort — missing fields simply become empty strings.
    """
    if not text:
        return None
    lines = text.splitlines()
    if not lines or not lines[0].startswith(">"):
        return None
    header = lines[0][1:]
    description = ""
    organism = ""
    # Description: everything after the 3rd '|' up to ' OS='
    parts = header.split("|", 3)
    tail = parts[3] if len(parts) >= 4 else header
    if " OS=" in tail:
        description = tail.split(" OS=", 1)[0].strip()
        rest = tail.split(" OS=", 1)[1]
        organism = rest.split(" OX=", 1)[0].strip() if " OX=" in rest else rest.strip()
    else:
        description = tail.strip()
    seq = "".join(line.strip() for line in lines[1:])
    if not seq:
        return None
    return seq, organism, description
