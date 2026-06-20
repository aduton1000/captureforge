# CaptureForge

**A reproducible Nextflow pipeline for hybridization-capture probe (bait) design — from a list of target genes to a synthesis-ready, QC'd bait set.**

![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
![Nextflow](https://img.shields.io/badge/Nextflow-DSL2%20%E2%89%A523.10-23aa62.svg)
![CLI: Python](https://img.shields.io/badge/CLI-Python-3776ab.svg)

## Overview

CaptureForge turns target genes and per-gene design intent into order-ready biotinylated-bait oligo sequences, a coverage BED, and a full QC report. It is disease- and organism-agnostic: point it at any gene list and a matching reference genome, and change inputs and config — never code. The scope is everything upstream of physical oligo synthesis; synthesis and wet-lab validation are out of scope.

CaptureForge is the **design** half of a pair. Its `final_covered_targets.bed` handoff feeds **[CallForge](https://github.com/aduton1000/callforge)**, the downstream germline-analysis pipeline.

## Features

- **End-to-end design** in ten Nextflow (DSL2) stages: target definition → masking → candidate tiling → specificity → composition → gap re-tiling → boost/dedup → synthesis scaffolding → QC aggregation → order package.
- **Design-class-aware**: every target interval carries a class (`coding`, `promoter`, `regulatory`, `cnv`, `str`, `anchor`) that drives tiling density, masking, GC policy, and boosting.
- **Specificity against the design genome**: mapping uses the no-alt primary assembly, enforced fail-loud so MHC and other alt-haplotype genes are not mistaken for off-target.
- **Reproducible**: pinned conda environment or containers, `provenance.json` + checksums on every run; references are fetched and checksummed explicitly, never downloaded mid-run.
- **Self-contained QC dashboard** plus per-gene coverage, GC/Tm, off-target, and copy-number figures.
- **The `captureforge` CLI** — `run`, `panel`, `index` — over a single Nextflow entrypoint, composed by profiles for laptop or cluster.

## Requirements

- **Nextflow** ≥ 23.10 (Java 17+). The repository bundles a launcher at `./.tools/nextflow`, so installing Nextflow separately is optional.
- **conda / mamba** for the pinned environment, **or** Docker / Apptainer for the container path.
- **Platforms**: Linux and macOS (Windows via WSL2). On Apple Silicon, build the environment under Rosetta with `CONDA_SUBDIR=osx-64`; the container path avoids this.

## Installation

```bash
git clone https://github.com/aduton1000/captureforge.git
cd captureforge

# pinned environment (creates the `captureforge` env)
mamba env create -f env/captureforge.yml                       # Linux
CONDA_SUBDIR=osx-64 mamba env create -f env/captureforge.yml   # macOS / Apple Silicon

# activate it so `captureforge` installs into this env, not your base environment
conda activate captureforge

# register the `captureforge` command (into the activated env)
pip install -e .
captureforge --version    # success: prints e.g. "captureforge 0.1.0"
```

A run composes **two** profiles with a comma: an **execution** profile — `local` (this machine: any Linux/macOS, Windows via WSL2) or `hpc_slurm` (a SLURM cluster) — and a **packaging** profile — `conda`, `docker`, or `apptainer` — that selects how dependencies are provided. For example, `-profile local,conda` runs on this machine using the conda environment.

## Quickstart

Run the bundled 10-gene test panel end-to-end. The fixture (a generic
`TP53/EGFR/BRCA1/…` mini-genome) ships pre-built in `test/data/`, so this works
straight from a clone — no build step, no local GRCh38 needed:

```bash
conda activate captureforge              # so the pipeline uses this env's tools
captureforge run -profile test           # or: -profile test,conda
bash test/check_outputs.sh               # sanity-check: expects all six design classes
```

This writes a complete design to `results/` — `baits.fasta` / `baits.csv`, `final_covered_targets.bed`, and `qc_dashboard.html` — in about two minutes. (The `test` profile uses a synthetic mini-genome fixture of generic, non-panel exemplar genes; it exercises the machinery, not real-panel quality. To rebuild/refresh the fixture for reproducibility, see `test/data/make_test_reference.sh` — not required to run the quickstart.)

## Usage

The three subcommands cover the full workflow; each is additive over the underlying Nextflow pipeline.

**Scaffold a panel spec** from a plain gene list into an editable CSV (`coding` is always implied; add `promoter`/`regulatory`/`cnv`/`str`/`anchor` intent per gene):

```bash
captureforge panel my_genes.txt -o my_panel.spec.csv
```

**Build the specificity index** once, from the no-alt primary assembly (the same genome the baits target). Use the **header-cleaned** FASTA that `bin/fetch_references.sh genome` writes (`*.clean.fa`, headers truncated to the first token) so the genome and the index agree on contig names:

```bash
captureforge index \
    --fasta refs/Homo_sapiens.GRCh38.dna.primary_assembly.clean.fa \
    --aligner blast \
    --out $PWD/refs/grch38_noalt
```

The same cleaned FASTA must be passed as `--genome_fasta`. (Stage 4 enforces the no-alt rule and a same-reference cross-check; it now also compares contig names on their first token, so an un-cleaned index still matches — but cleaning at prep time is the supported path.)

**Run the design** on a real panel:

```bash
captureforge run -profile local,conda \
    --input my_panel.spec.csv \
    --genome_fasta refs/GRCh38.noalt.fa \
    --gtf refs/GRCh38.112.gtf.gz \
    --specificity_db $PWD/refs/grch38_noalt \
    --outdir results
```

`--specificity_db` must be an absolute path to the index prefix, built from the no-alt primary assembly; the pipeline fails loud on a with-alt or toplevel index.

**On a shared SLURM cluster**, install once centrally and run through the `captureforge-run` wrapper. A site fills in `share/captureforge-env.sh.example` once (PATH, image, refs, container bind) and sources it from `~/.bashrc`, and copies `conf/cluster.params.example.yaml` once; then every run is just:

```bash
captureforge-run -params-file site.params.yaml --input panel.txt --outdir results
```

The container runs under **Singularity or Apptainer** (same `.sif`); references that live outside `$HOME`/CWD (the specificity index) are made visible by binding their root (`CAPTUREFORGE_BIND`), and bowtie2 specificity uses capped `-k` reporting so it scales to the whole genome. The full step-by-step — build, reference prep (incl. header cleaning), the env file, and the bind — is in [`docs/hpc_deployment.md`](docs/hpc_deployment.md) and the manual.

## Inputs and outputs

**Inputs.** A panel definition — a plain gene list, a CSV/TSV spec (`gene` + optional `promoter,regulatory,cnv,str,extra_regions,anchors`), or a curated Excel workbook — plus a no-alt reference FASTA, a GTF, and a specificity index built from that FASTA. References are managed by `bin/fetch_references.sh` (fetch + checksum + provenance), never fetched mid-run. Required parameters: `--input`, `--genome_fasta`, `--gtf`.

**Outputs (`results/`).** `baits.fasta` and `baits.csv` (order-ready sequences + metadata); `final_covered_targets.bed` (the handoff to CallForge); `qc_dashboard.html` (self-contained QC report); `metrics.json` and per-stage TSVs; `MANIFEST` + `provenance.json` (checksums, tool/reference/parameter provenance).

## Documentation

The full user manual covers concepts, laptop and shared-HPC installs, the start-from-a-gene-name tutorial, running on other organisms, QC interpretation, and troubleshooting:

- **Manual**: [`docs/manual/captureforge_manual.md`](docs/manual/captureforge_manual.md) · [PDF](docs/manual/captureforge_manual.pdf)
- **HPC deployment**: [`docs/hpc_deployment.md`](docs/hpc_deployment.md)
- **Outputs overview**: [`docs/outputs_overview.md`](docs/outputs_overview.md) — what a run produces (bait set, the as-built `final_covered_targets.bed`, and the QC dashboard). The manual works a full example end-to-end starting from a single gene name, and ships a `test`-profile run against a synthetic mini-genome fixture.

## CaptureForge and CallForge

CaptureForge (design) and **[CallForge](https://github.com/aduton1000/callforge)** (analysis) are a pair: CaptureForge defines, tiles, and QCs the capture baits and emits `final_covered_targets.bed` plus per-gene metadata; CallForge consumes that handoff to align, call, annotate, and run cohort QC and rare-variant burden testing. Each runs independently.

## Citation, license, and issues

- **Cite**: see [`CITATION.cff`](CITATION.cff).
- **License**: MIT — see [`LICENSE`](LICENSE).
- **Issues / contributions**: open an issue or pull request on the [GitHub repository](https://github.com/aduton1000/captureforge).
