# CaptureForge — gene-list input schema

CaptureForge auto-detects the input format from its extension/content and normalizes **any** of the three
on-ramps below into a single internal **target table** that the rest of the pipeline consumes. No prior
project file or script is ever required.

```
   plain text  ─┐
   CSV/TSV spec ─┼─▶  bin/parse_gene_list.py  ─▶  targets.tsv (internal table)  ─▶  Stage 1…10
   Excel book  ─┘
```

## On-ramp 1 — Plain-text symbol list  (`.txt`, `.list`, or one column)

One gene symbol per line. Blank lines and `#` comments ignored. Applies the **default design classes**
(`params.default_classes`, default `coding,promoter`). Zero friction — the way to run with *nothing but a
list of genes*.

```
GENE
GENE
GENE
```

## On-ramp 2 — CSV/TSV spec  (`.csv`, `.tsv`)

Header **`gene`** required; all others optional. Boolean-ish columns accept `1/0`, `yes/no`, `true/false`.

| column | type | meaning |
|---|---|---|
| `gene` | string (**required**) | gene symbol (aliases in parentheses tolerated) |
| `promoter` | bool | add strand-aware promoter window |
| `regulatory` | bool | intersect Ensembl Regulatory / cCREs for distal elements |
| `cnv` | bool | add even-tiled CNV windows |
| `str` | bool | tile unique flanks only (skip repeat core) |
| `extra_regions` | string | extra intervals `chr:start-end` separated by `;` |
| `anchors` | string | rsIDs / variant tokens separated by `,` or `;` |

`coding` is always implied. Example: `assets/gene_lists/example_spec.csv`.

## On-ramp 3 — Rich curated workbook  (`.xlsx`)

Reads the sheet named by `--excel-sheet` (default `Candidate Panel (expanded)`). The header row is detected
automatically (the row containing a `Gene` cell). Recognized columns:

| workbook column | used for |
|---|---|
| `Gene` (**required**) | symbol (alias-in-parentheses tolerated, e.g. `GENE (GENE)`) |
| `Variant class` | mapped to design classes (see mapping below) |
| `Key variants to capture` | scanned for `rs\d+` → anchor targets |
| `Locus` | cytoband (informational / QC) |
| `Burden group` | rare-variant set-test grouping (carried through, informational) |

### `Variant class` → design-class mapping (case-insensitive substring match)

| text contains | adds class |
|---|---|
| `reg`, `promoter`, `enhancer`, `regulatory`, `haplotype` | `promoter` (+ `regulatory` if `enhancer`/`distal`/`regulatory`) |
| `cnv`, `copy number`, `structural`, `deletion`, `duplication` | `cnv` |
| `str`, `repeat`, `tandem`, `microsatellite`, `expansion` | `str` |
| `coding`, `LoF`, `burden`, `missense`, `indel`, anything else | `coding` (always added) |

## Normalization rules (all on-ramps)

1. **Alias resolution.** A symbol written `A (B)` is resolved to whichever of `A`/`B` is a valid current
   symbol in the supplied GTF (`gene_name`). Examples: `GENE (GENE)`→`GENE`, `GENE (GENE)`→`GENE`,
   `GENE (GENE)`→`GENE`, `GENE (GENE)`→`GENE`. Case-insensitive fallback is attempted.
2. **Validation.** Every symbol must resolve to an Ensembl gene in the GTF. Unresolvable symbols cause a
   **loud failure** with the offending list (override with `--allow-unresolved` to skip+warn).
3. **Anchors.** `rs\d+` tokens are extracted as anchor targets (±`anchor_flank` bp, default 150).

## Internal target table (`targets.tsv`) — output of `parse_gene_list.py`

Tab-separated, one row per gene; everything downstream consumes only this file:

| column | example |
|---|---|
| `gene` | `GENE` |
| `ensembl_id` | `ENSG00000141510` |
| `chrom` | `1` |
| `strand` | `-` |
| `gene_start` | `7668421` |
| `gene_end` | `7687550` |
| `design_classes` | `coding|promoter|regulatory` |
| `anchors` | `rs1042522,rs1042522` |
| `extra_regions` | `1:203650000-203651000` |
| `burden_group` | `single` |
