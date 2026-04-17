"""gnomAD constraint client.

gnomAD (Genome Aggregation Database) provides population genetics data from
>125,000 exomes and >15,000 whole genomes.  This client fetches gene-level
constraint metrics that quantify how much a gene tolerates loss-of-function
(LoF) and missense mutations in the human population — a critical pre-filter
for protein engineering campaigns.

Key metrics:
    pLI     — probability of being loss-of-function intolerant (>0.9 = intolerant)
    oe_lof  — observed/expected LoF variant ratio (lower = more constrained)
    LOEUF   — oe_lof_upper confidence interval bound (most used constraint metric)
    oe_mis  — observed/expected missense ratio
    lof_z   — LoF Z-score (>3 = significant constraint)
    mis_z   — missense Z-score

Data source: gnomAD v4 via public GraphQL API.  No API key required.
"""

from __future__ import annotations

import asyncio
import logging

import httpx

from genesis_bio_mcp.models import GnomADConstraint

logger = logging.getLogger(__name__)

_GNOMAD_API = "https://gnomad.broadinstitute.org/api"
_SEMAPHORE = asyncio.Semaphore(3)

_CONSTRAINT_QUERY = """
query GeneConstraint($symbol: String!) {
  gene(gene_symbol: $symbol, reference_genome: GRCh38) {
    gene_id
    name
    canonical_transcript_id
    gnomad_constraint {
      pLI
      lof_z
      mis_z
      oe_lof
      oe_lof_lower
      oe_lof_upper
      oe_mis
      exp_lof
      exp_mis
      obs_lof
      obs_mis
    }
  }
}
"""

# Fetch all variants for a gene with enough metadata to pick the matching
# protein change. gnomAD's `gene.variants` returns every variant catalogued
# against the canonical transcript.
_GENE_VARIANTS_QUERY = """
query GeneVariants($symbol: String!) {
  gene(gene_symbol: $symbol, reference_genome: GRCh38) {
    variants(dataset: gnomad_r4) {
      variant_id
      hgvsp
      consequence
    }
  }
}
"""


class GnomADClient:
    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client
        self._cache: dict[str, GnomADConstraint] = {}
        # Per-gene variants cache. The payload is small (tens of thousands of
        # IDs for a large gene) so one fetch per gene per session is fine.
        self._variants_cache: dict[str, list[dict] | None] = {}

    async def get_constraint(self, gene_symbol: str) -> GnomADConstraint | None:
        """Return gnomAD constraint metrics for *gene_symbol*.

        Returns ``None`` if the gene is not found or the API is unavailable.
        Returns a :class:`GnomADConstraint` with ``constraint_available=False``
        if the gene exists but has no constraint data (e.g. insufficiently covered).
        """
        symbol = gene_symbol.strip().upper()
        if symbol in self._cache:
            logger.debug("gnomAD cache hit: %s", symbol)
            return self._cache[symbol]

        async with _SEMAPHORE:
            result = await self._fetch(symbol)

        if result is not None:
            self._cache[symbol] = result
        return result

    async def _fetch(self, symbol: str) -> GnomADConstraint | None:
        try:
            resp = await self._client.post(
                _GNOMAD_API,
                json={"query": _CONSTRAINT_QUERY, "variables": {"symbol": symbol}},
                timeout=20.0,
            )
            resp.raise_for_status()
            data = resp.json()
        except Exception as exc:
            logger.warning("gnomAD fetch failed for %s: %s", symbol, exc)
            return None

        errors = data.get("errors")
        if errors:
            logger.warning("gnomAD GraphQL error for %s: %s", symbol, errors)
            return None

        gene = (data.get("data") or {}).get("gene")
        if not gene:
            logger.info("gnomAD: no gene record for '%s'", symbol)
            return None

        constraint = gene.get("gnomad_constraint")
        if not constraint:
            return GnomADConstraint(
                gene_symbol=symbol,
                ensembl_id=gene.get("gene_id"),
                gene_name=gene.get("name"),
                constraint_available=False,
            )

        def _f(key: str) -> float | None:
            v = constraint.get(key)
            return float(v) if v is not None else None

        def _i(key: str) -> int | None:
            v = constraint.get(key)
            return int(v) if v is not None else None

        return GnomADConstraint(
            gene_symbol=symbol,
            ensembl_id=gene.get("gene_id"),
            gene_name=gene.get("name"),
            constraint_available=True,
            pLI=_f("pLI"),
            lof_z=_f("lof_z"),
            mis_z=_f("mis_z"),
            oe_lof=_f("oe_lof"),
            oe_lof_lower=_f("oe_lof_lower"),
            oe_lof_upper=_f("oe_lof_upper"),
            oe_mis=_f("oe_mis"),
            exp_lof=_f("exp_lof"),
            exp_mis=_f("exp_mis"),
            obs_lof=_i("obs_lof"),
            obs_mis=_i("obs_mis"),
        )

    async def find_variant_id_by_protein_change(
        self, gene_symbol: str, hgvs_protein: str
    ) -> str | None:
        """Look up the gnomAD variant_id (chrom-pos-ref-alt) for a protein change.

        gnomAD does not expose a direct protein → variant endpoint; we load
        the gene's full variant list once per session and filter on the
        canonical ``hgvsp`` string.

        Args:
            gene_symbol: HGNC symbol, e.g. ``"TP53"``.
            hgvs_protein: Canonical HGVS p. notation, e.g. ``"p.Arg175His"``.

        Returns:
            The variant_id string on success, or ``None`` if the gene /
            variant is not indexed or the API is unreachable.
        """
        symbol = gene_symbol.strip().upper()
        target = hgvs_protein.strip()
        variants = await self._load_variants(symbol)
        if not variants:
            return None
        for v in variants:
            hgvsp = v.get("hgvsp")
            if hgvsp and hgvsp.strip() == target:
                return v.get("variant_id")
        return None

    async def _load_variants(self, symbol: str) -> list[dict] | None:
        if symbol in self._variants_cache:
            logger.debug("gnomAD variants cache hit: %s", symbol)
            return self._variants_cache[symbol]
        async with _SEMAPHORE:
            try:
                resp = await self._client.post(
                    _GNOMAD_API,
                    json={"query": _GENE_VARIANTS_QUERY, "variables": {"symbol": symbol}},
                    timeout=30.0,
                )
                resp.raise_for_status()
                data = resp.json()
            except Exception as exc:
                logger.warning("gnomAD variants fetch failed for %s: %s", symbol, exc)
                self._variants_cache[symbol] = None
                return None
        gene = (data.get("data") or {}).get("gene")
        variants = (gene or {}).get("variants") if gene else None
        if not isinstance(variants, list):
            self._variants_cache[symbol] = None
            return None
        self._variants_cache[symbol] = variants
        return variants
