"""UniProt REST API client for protein annotation retrieval."""

from __future__ import annotations

import logging
from typing import Optional

import httpx

from genesis_bio_mcp.models import KnownVariant, ProteinInfo

logger = logging.getLogger(__name__)

_BASE_URL = "https://rest.uniprot.org/uniprotkb"
_FIELDS = (
    "accession,gene_names,protein_name,organism_name,"
    "cc_function,cc_subcellular_location,cc_disease,"
    "xref_reactome,xref_pdb,ft_variant"
)


class UniProtClient:
    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client

    async def get_protein(self, gene_symbol: str) -> Optional[ProteinInfo]:
        """Return Swiss-Prot annotation for a human gene symbol, or None if not found."""
        symbol = gene_symbol.strip().upper()
        data = await self._search(symbol, reviewed_only=True)
        if data is None:
            data = await self._search(symbol, reviewed_only=False)
        if data is None:
            return None
        return _parse_entry(data, symbol)

    async def search_by_synonym(self, synonym: str) -> Optional[dict]:
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

    async def _search(self, symbol: str, *, reviewed_only: bool) -> Optional[dict]:
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
            logger.warning("UniProt search failed for '%s' (reviewed=%s): %s", symbol, reviewed_only, exc)
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

    # Variants
    features = entry.get("features", [])
    known_variants: list[KnownVariant] = []
    for feat in features:
        if feat.get("type") != "Natural variant":
            continue
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
        reviewed=reviewed,
    )
