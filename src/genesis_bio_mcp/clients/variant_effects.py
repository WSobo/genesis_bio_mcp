"""Aggregator for mutation-level variant effect lookups.

Composes existing clients (gnomAD for variant-ID resolution, MyVariant.info
for ClinVar + AlphaMissense + REVEL + CADD + gnomAD frequency, MaveDB for
per-variant DMS fitness scores) into a single :class:`VariantEffects`
report.

Kept lightweight — no HTTP of its own; all network traffic flows through
the underlying clients (each of which has its own semaphore, cache, and
error handling).
"""

from __future__ import annotations

import asyncio
import logging

from genesis_bio_mcp.clients.gnomad import GnomADClient
from genesis_bio_mcp.clients.mavedb import MaveDBClient
from genesis_bio_mcp.clients.myvariant import MyVariantClient
from genesis_bio_mcp.models import MaveDBVariantScore, VariantEffects
from genesis_bio_mcp.tools.variant_parser import (
    canonical_one_letter,
    canonical_three_letter,
    parse_protein_change,
)

logger = logging.getLogger(__name__)

# How many DMS score sets to probe per-variant. The top-N by variant count
# cover the high-coverage datasets; lower-coverage sets rarely include a
# specific variant and cost one HTTP call each to check.
#
# Tuning note: raising this to 5 roughly doubles latency for genes with
# many DMS datasets (e.g. BRCA1 has 30+), for marginal recall gain since
# the largest sets are saturating. Lowering below 2 risks missing fitness
# data for genes where the largest MaveDB entry doesn't cover the variant.
_MAX_DMS_SCORESETS_TO_PROBE = 3


class VariantEffectsClient:
    """Aggregator client — no state of its own beyond injected clients."""

    def __init__(
        self,
        *,
        gnomad: GnomADClient,
        myvariant: MyVariantClient,
        mavedb: MaveDBClient,
    ) -> None:
        self._gnomad = gnomad
        self._myvariant = myvariant
        self._mavedb = mavedb

    async def get_effects(self, gene_symbol: str, mutation: str) -> VariantEffects:
        """Return a consolidated variant-effect report for *gene_symbol* + *mutation*.

        Args:
            gene_symbol: HGNC symbol. Case-insensitive.
            mutation: Protein change string. See
                :func:`genesis_bio_mcp.tools.variant_parser.parse_protein_change`
                for accepted forms.

        Returns:
            :class:`VariantEffects`. Fields that could not be resolved are
            ``None`` or empty lists with an explanatory entry in ``notes``.
        """
        symbol = gene_symbol.strip().upper()
        orig, pos, new = parse_protein_change(mutation)
        one_letter = canonical_one_letter(orig, pos, new)
        hgvs_p = canonical_three_letter(orig, pos, new)

        notes: list[str] = []

        # Step 1: resolve variant_id via gnomAD (single HTTP call that loads
        # all the gene's variants, then in-memory filter). Somatic cancer
        # hotspots (BRAF V600E, KRAS G12C/D, EGFR L858R, ...) are *expected*
        # to be absent from gnomAD — the population reference is germline
        # only. Downstream ClinVar/AlphaMissense lookups must NOT be gated
        # on this; we fall back to a protein-change query in step 2.
        variant_id = await self._gnomad.find_variant_id_by_protein_change(symbol, hgvs_p)
        if variant_id is None:
            notes.append(
                f"{symbol} {one_letter} was not found in gnomAD v4 — population frequency "
                "unavailable (typical for somatic cancer variants)."
            )

        # Step 2: fan out — MyVariant.info (annotation) and MaveDB (DMS) in
        # parallel. When gnomAD has no variant_id, fall back to MyVariant's
        # /query endpoint keyed on gene + protein position + alt AA so
        # ClinVar/AlphaMissense/CADD are still retrieved.
        if variant_id:
            annotation_task = self._myvariant.get_annotation(
                _variant_id_to_hgvs_genomic(variant_id)
            )
        else:
            _, pos, new = parse_protein_change(mutation)
            annotation_task = self._myvariant.query_by_protein_change(symbol, pos, new)
        dms_task = self._collect_dms_scores(symbol, hgvs_p)

        annotation, dms_scores = await asyncio.gather(annotation_task, dms_task)

        if annotation is None:
            notes.append(
                "MyVariant.info returned no record for this variant — it may be too rare to be "
                "indexed with downstream annotations yet."
            )

        return VariantEffects(
            gene_symbol=symbol,
            mutation_input=mutation,
            canonical_one_letter=one_letter,
            canonical_hgvs_protein=hgvs_p,
            gnomad_variant_id=variant_id,
            annotation=annotation,
            dms_scores=dms_scores,
            notes=notes,
        )

    async def _collect_dms_scores(self, gene_symbol: str, hgvs_p: str) -> list[MaveDBVariantScore]:
        """Probe top DMS score sets for a per-variant fitness score."""
        results = await self._mavedb.get_dms_scores(gene_symbol)
        if results is None or not results.score_sets:
            return []
        top = results.score_sets[:_MAX_DMS_SCORESETS_TO_PROBE]
        per_set = await asyncio.gather(
            *[self._mavedb.get_variant_score(ss.urn, hgvs_p, ss.title) for ss in top],
            return_exceptions=True,
        )
        flattened: list[MaveDBVariantScore] = []
        for entry in per_set:
            if isinstance(entry, Exception):
                logger.warning("MaveDB per-variant probe failed: %s", entry)
                continue
            flattened.extend(entry)
        return flattened


def _variant_id_to_hgvs_genomic(variant_id: str) -> str:
    """Translate gnomAD variant_id to MyVariant.info's expected HGVS form.

    gnomAD: ``"17-7675088-C-T"`` → MyVariant: ``"chr17:g.7675088C>T"``.
    """
    chrom, pos, ref, alt = variant_id.split("-")
    return f"chr{chrom}:g.{pos}{ref}>{alt}"
