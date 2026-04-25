# Protein engineering workflows

genesis-bio-mcp v0.2.0 added three protein-engineering tools that compose
into realistic design workflows: sequence-level analysis, mutation-level
pathogenicity, and T-cell immunogenicity.

This document walks through the three core workflows with concrete inputs
and expected outputs. See [tools.md](tools.md) for the full tool reference
and [architecture.md](architecture.md) for client design.

- [Sequence-level analysis](#sequence-level-analysis)
- [Mutation-level pathogenicity and fitness](#mutation-level-pathogenicity-and-fitness)
- [Therapeutic immunogenicity / T-cell epitope screening](#therapeutic-immunogenicity--t-cell-epitope-screening)
- [Combined workflow example](#combined-workflow-example)

---

## Sequence-level analysis

**Tool:** `get_protein_sequence`

**When to call:** before any mutagenesis, region-selective design, antibody
CDR work, or developability review. Gets the FASTA and computes ExPASy
ProtParam-equivalent biochemistry plus a liability-motif scan in one call.

### What it returns

- UniProt accession + organism + protein name
- Full sequence (optionally sliced to a residue range)
- **Biochem features**
  - Molecular weight (Da)
  - Theoretical pI (Bjellqvist pKa set — matches ExPASy ProtParam to ±0.1 pH)
  - GRAVY (Kyte-Doolittle mean hydropathy)
  - Net charge at pH 7.4
  - Aromatic fraction (F + W + Y)
  - Cysteine count
  - ε₂₈₀ reduced and oxidized (Edelhoch / Pace formula)
- **Liability motif hits**, grouped by class:
  - `deamidation`  — NG, NS (Asn → Asp/isoAsp shelf-life degradation)
  - `isomerization`  — DG, DS (Asp → isoAsp via succinimide)
  - `n_glycosylation`  — `N[^P][ST]` sequon (glycan attachment site)
  - `oxidation_methionine`, `oxidation_tryptophan`  — positions as data, not flagged risks
  - `free_cysteine`  — only emitted when UniProt DISULFID annotation confirms the Cys is **not** bonded (never on parity heuristic alone)
  - `cysteine_position`  — Cys positions reported as data when no DISULFID annotation exists
- UniProt-annotated disulfide bond positions (if any)

### Example

```text
get_protein_sequence(gene_symbol="BRAF", start=450, end=620)

→ Length: 171 aa
→ MW: 19,358.4 Da | pI: 9.24 | GRAVY: −0.155 | net charge @ pH 7.4: +3.83
→ ε₂₈₀ reduced: 33,460 M⁻¹cm⁻¹ (no disulfides in this region)
→ Liability motifs
  - isomerization (1): [5]
  - n_glycosylation (1): [37]
  - oxidation_methionine (5): [35, 68, 101, 115, 171]
  - oxidation_tryptophan (5): [1, 27, 82, 155, 170]
  - cysteine_position (1): [83]
```

Documentation for every constant (AA masses, pKa values, hydropathy scale,
extinction coefficients) is inline in
`src/genesis_bio_mcp/tools/biochem.py` with source citations.

---

## Mutation-level pathogenicity and fitness

**Tool:** `get_variant_effects`

**When to call:** to assess whether a specific mutation (engineered or
observed) is likely pathogenic, tolerated, or functionally disruptive.
Fans out to three sources in parallel.

### Pipeline

```
       parse R175H → p.Arg175His
              │
              ▼
  ┌──────────────────────────┐
  │ gnomAD v4 (GraphQL)      │  resolve protein change → variant_id
  │ find_variant_id_by_hgvsp │  e.g. "17-7675088-C-T"
  └──────────┬───────────────┘
             │
             ▼
  ┌─────────────────────────────────────────┐   ┌───────────────────────┐
  │ MyVariant.info (hg38 genomic HGVS key)  │   │ MaveDB (per-variant)  │
  │                                         │   │                       │
  │ • ClinVar RCV records                   │   │ probes top-3 DMS      │
  │ • AlphaMissense score + class           │   │ score sets for the    │
  │ • REVEL / CADD / SIFT / PolyPhen-2      │   │ gene, filters by      │
  │ • gnomAD exome AF + population split    │   │ hgvs_pro match        │
  └─────────────────────────────────────────┘   └───────────────────────┘
             │                                           │
             └───────────── asyncio.gather ──────────────┘
                                 │
                                 ▼
                         VariantEffects report
```

### Example

```text
get_variant_effects(gene_symbol="TP53", mutation="R175H")

→ TP53 R175H → p.Arg175His → gnomAD variant 17-7675088-C-T

ClinVar: Pathogenic (14 submissions) | AlphaMissense: 0.97 (likely_pathogenic)
DMS: 2 per-variant scores across 2 score sets

gnomAD exome AF: 3.98×10⁻⁶ (af_nfe=8.80×10⁻⁶, other populations = 0)

In silico predictions
| Predictor    | Score  | Interpretation       |
|--------------|--------|----------------------|
| AlphaMissense| 0.973  | likely_pathogenic    |
| REVEL        | 0.922  | likely pathogenic    |
| CADD (Phred) | 25.9   | top 1% deleterious   |
| SIFT         | 0.000  | deleterious          |

ClinVar submissions (top 5)
- Pathogenic (RCV000013173) — Li-Fraumeni syndrome 1 | Review: 4 submitters, no conflicts
- Pathogenic (RCV000131301) — Hereditary cancer-predisposing syndrome
- Pathogenic (RCV000204931) — Li-Fraumeni syndrome | Review: expert panel
- ...and 11 more

MaveDB DMS scores
- urn:mavedb:00000068-a-1: 1.79   — Mutated p53 paired with wildtype under nutlin-3
- urn:mavedb:00000068-c-1: −0.74  — Mutated p53 with etoposide
```

### Accepted mutation formats

`parse_protein_change` normalizes: `R175H`, `p.R175H`, `Arg175His`,
`p.Arg175His` (case-insensitive). Anything else raises `ValueError` with a
message listing the accepted forms — invalid input surfaces as a tool
error, never silent coercion.

### Splice / UTR / regulatory consequences

For non-missense variants — splice-site, 5'/3' UTR, intronic, regulatory —
use `get_variant_consequences` (Ensembl VEP). It complements
`get_variant_effects` (which sources dbNSFP missense scores via
MyVariant.info) with VEP's transcript-aware consequence terms and its own
SIFT/PolyPhen calls. Canonical transcript by default;
`include_all_transcripts=True` for the full set. Inside
`get_variant_effects`, VEP is already fanned out as a third parallel
source, so a separate call is only needed when you want the
all-transcripts view or are starting from coordinates rather than a
protein change.

---

## Therapeutic immunogenicity / T-cell epitope screening

**Tool:** `get_mhc_binding`

**When to call:** for any therapeutic protein (antibody CDR, enzyme,
vaccine antigen) to predict which peptides would be presented by MHC and
could trigger T-cell recognition.

### Defaults

| Parameter | Default | Notes |
|---|---|---|
| `mhc_class` | `"I"` | `"II"` for CD4 / helper-T presentation |
| `hla_alleles` | 5-allele IEDB class-I reference | `HLA-A*02:01`, `A*03:01`, `A*24:02`, `B*07:02`, `B*35:01` — ~85% global phenotype coverage |
| Class II default panel | 5-allele DRB1 | `DRB1*01:01`, `*03:01`, `*04:01`, `*07:01`, `*15:01` — ~60% coverage |
| `peptide_lengths` | `[9, 10]` for class I, `[15]` for class II | Canonical per IEDB |
| `method` | `netmhcpan_el` (class I), `netmhciipan_el` (class II) | IEDB-recommended since 2023.09 |

### IEDB binder-class thresholds

- Strong binder: percentile rank **< 0.5** — high-confidence binding
- Weak binder: **< 2.0** — likely binding
- Non-binder: **≥ 2.0**

### Implementation notes

- The NextGen Tools service uses an async ticket pattern: the client POSTs
  the pipeline spec, polls the result URI every 2 s up to a 60 s deadline,
  and parses the completed `peptide_table` block.
- Whole proteins are auto-windowed by IEDB into peptides of the requested
  length; a `sequence` longer than `peptide_lengths[0]` will generate
  multiple windows.
- Request-size cap (peptides × alleles ≤ 2000) raises `ValueError` early
  so accidental whole-proteome submissions fail fast.

### Example — HIV Gag SLYNTVATL (known A\*02:01-restricted CTL epitope)

```text
get_mhc_binding(
  sequence="SLYNTVATL",
  alleles=["HLA-A*02:01", "HLA-A*03:01", "HLA-B*07:02"],
  mhc_class="I",
)

→ 1 strong binder, 1 weak, across 3 peptide×allele rows

| Peptide    | Allele       | %tile | Score | Class       |
|------------|--------------|-------|-------|-------------|
| SLYNTVATL  | HLA-A*02:01  | 0.06  | 0.828 | strong      |
| SLYNTVATL  | HLA-B*07:02  | 1.80  | 0.045 | weak        |
| SLYNTVATL  | HLA-A*03:01  | 3.00  | 0.028 | non_binder  |

Strong binders per allele: HLA-A*02:01=1
```

The A\*02:01 strong-binder classification matches the canonical
immunodominance literature for this epitope.

---

## Combined workflow example

A realistic therapeutic-engineering check combining all three tools.
All of these can also be driven by a free-text prompt to
`run_biology_workflow`, which chains them automatically.

### "Flag any liabilities in the trastuzumab HC CDR3 and check its T-cell immunogenicity risk."

1. `get_protein_sequence` — scan the CDR3 sequence for deamidation /
   isomerization / oxidation hotspots and get charge + pI for developability.
2. `get_mhc_binding(sequence=CDR3, mhc_class="I")` — check that no 9mer
   window in CDR3 is a strong binder to common HLA alleles.
3. If the antibody targets a human antigen, also:
   `get_variant_effects` on any paratope-contact residue variant that
   might affect binding.

### "Assess the engineered mutation BRCA1 C47R for clinical risk and DMS fitness."

1. `get_variant_effects(gene_symbol="BRCA1", mutation="C47R")` returns the
   ClinVar verdict, AlphaMissense pathogenicity, REVEL/CADD/SIFT/PolyPhen,
   gnomAD frequency, and any matching MaveDB DMS fitness score — one call.
2. `get_variant_constraints("BRCA1")` for gene-level LoF / missense
   tolerance (pLI, LOEUF) — sanity check that BRCA1 is constrained (it is).
3. `get_dms_scores("BRCA1")` lists every available DMS score set (>30 for
   BRCA1) if deeper saturation data is needed.

### "Find potential T-cell epitopes in a therapeutic enzyme before humanization."

1. `get_protein_sequence` to retrieve the full sequence and flag
   N-glycosylation sites that can mask epitopes.
2. `get_mhc_binding(sequence=full_protein, mhc_class="II")` to screen
   15mer windows against the default DRB1 panel. Strong binders
   (percentile < 0.5) are the highest-priority targets for humanizing
   mutations.
3. Cross-reference with `get_epitope_data` for known B-cell epitopes on
   the same protein so you understand both arms of the adaptive response.

---

## What's next (deferred)

v0.2.0 shipped M1–M3 of the protein-engineering roadmap. Future batches
(tracked in the project plan):

- **M4** — enhance `get_antibody_structures` with CDR developability
  metrics (charge, pI, GRAVY, liability motifs on CDRs) — reuses `biochem.py`
- **M5** — `get_structural_homologs` via Foldseek web API
- **M6** — per-residue pLDDT bands / disordered-region detection in
  `get_protein_structure` from AlphaFold confidence JSON
- **M7** — promote the MaveDB per-variant extension in `get_variant_effects`
  to a first-class `get_dms_variant_score` tool

Explicitly out of scope (APIs-only constraint): ESMFold, ProteinMPNN,
Rosetta / FoldX, ANARCI / AbNumber, MHCflurry — all require local GPU or
heavyweight installs.
