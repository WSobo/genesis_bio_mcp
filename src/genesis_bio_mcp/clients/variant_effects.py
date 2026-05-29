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

from genesis_bio_mcp.clients.ensembl import EnsemblClient
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


class VariantEffectsClient:
    """Aggregator client — no state of its own beyond injected clients."""

    def __init__(
        self,
        *,
        gnomad: GnomADClient,
        myvariant: MyVariantClient,
        mavedb: MaveDBClient,
        ensembl: EnsemblClient,
    ) -> None:
        self._gnomad = gnomad
        self._myvariant = myvariant
        self._mavedb = mavedb
        self._ensembl = ensembl

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
        # VEP consequences are orthogonal to MyVariant's dbNSFP-sourced SIFT/
        # PolyPhen — Ensembl returns its own predictors plus splice/UTR/
        # regulatory overlap, which dbNSFP doesn't cover.
        vep_task = self._ensembl.get_vep_consequences(symbol, hgvs_p)

        annotation, dms_scores, vep = await asyncio.gather(annotation_task, dms_task, vep_task)

        if annotation is None:
            notes.append(
                "MyVariant.info returned no record for this variant — it may be too rare to be "
                "indexed with downstream annotations yet."
            )
        if vep is None:
            notes.append(
                "Ensembl VEP returned no consequence prediction — canonical transcript "
                "may be unavailable or the HGVS form was rejected."
            )

        return VariantEffects(
            gene_symbol=symbol,
            mutation_input=mutation,
            canonical_one_letter=one_letter,
            canonical_hgvs_protein=hgvs_p,
            gnomad_variant_id=variant_id,
            annotation=annotation,
            dms_scores=dms_scores,
            vep_consequences=vep,
            notes=notes,
        )

    async def _collect_dms_scores(self, gene_symbol: str, hgvs_p: str) -> list[MaveDBVariantScore]:
        """Probe top DMS score sets for a per-variant fitness score.

        Delegates to :meth:`MaveDBClient.get_variant_scores_for_gene` so the
        gene→score-sets→per-variant probe logic lives in one place (also used by
        the standalone ``get_dms_variant_score`` tool).
        """
        return await self._mavedb.get_variant_scores_for_gene(gene_symbol, hgvs_p)


def _variant_id_to_hgvs_genomic(variant_id: str) -> str:
    """Translate gnomAD variant_id to MyVariant.info's expected HGVS form.

    gnomAD: ``"17-7675088-C-T"`` → MyVariant: ``"chr17:g.7675088C>T"``.
    """
    chrom, pos, ref, alt = variant_id.split("-")
    return f"chr{chrom}:g.{pos}{ref}>{alt}"
