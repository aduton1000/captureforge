# What CaptureForge produces

This page summarizes the deliverables of a CaptureForge run at the method level. It describes the
**shape** of the outputs, not any specific panel — for how to read and act on them, see the
[user manual](manual/captureforge_manual.md) (chapters *Outputs explained* and *Interpreting QC and
coverage*).

Every run writes to `results/<outdir>/`. Three deliverables matter to most users.

## 1. A synthesis-ready bait set — `order_package/`

- **`baits.fasta`** — the final, order-ready bait constructs. For the default `ivt_rna` format each oligo is
  `handle5 + T7 + bait + handle3`; `dna_pool` and `ready_to_use` emit the bait core only.
- **`baits.csv`** — one row per unique bait with full metadata: stable `bait_id`, sequence, copy number,
  source gene(s) and design class(es), source region, strand, GC, Tm, and an off-target flag.
- **`order_generic.tsv`** — a vendor-neutral order table (`Pool_Name, Probe_ID, Sequence, Replicates`) to map
  onto a synthesizer's template.
- **`MANIFEST`** + **`provenance.json`** — md5 checksums plus the tool versions, reference releases, git
  commit, design genome, and specificity reference that produced the run.

## 2. The as-built coverage handoff — `stage9_qc/final_covered_targets.bed`

The **intersection of the raw targets with what the final baits actually cover**, labelled `gene|class`. This
is the AS-BUILT footprint — the regions you can expect reads in — and is the handoff to a downstream
variant-analysis pipeline. Measure downstream calls against this file, not the nominal target set. Alongside
it, `metrics.json`, `coverage_by_gene.tsv`, and `coverage_by_class.tsv` carry the machine-readable QC.

## 3. A self-contained QC dashboard — `qc_dashboard.html` + `figures/`

A single HTML file (open in a browser) embeds every figure and the key metric tables: coverage by class and
by gene, GC/Tm composition, the off-target distribution, copy-number boosting, re-tiling convergence, the
cost-vs-tier check, and genome-browser-style per-gene track panels. All figures are also written as PNG + SVG
under `figures/`.

Per-stage intermediates live under `stage1_targets/` … `stage8_synthesis/`, and `pipeline_info/` holds the
Nextflow `report` / `timeline` / `trace` / `dag`. The bundled `test` profile produces this full set against a
synthetic mini-genome fixture, which is the basis for the illustrative figures in the manual.
