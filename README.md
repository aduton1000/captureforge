# CaptureForge

**A general-purpose, disease- and organism-agnostic hybridization-capture probe (bait) design pipeline.**

CaptureForge turns a **gene list + per-gene design intent** into **order-ready biotinylated-bait
oligo sequences** plus a full QC/design report. It is a reusable *method tool*: the first use is a
multi-gene human study-susceptibility panel, but nothing about the structure or defaults is
study-specific — point it at any gene list (any disease, any study) and, where feasible, any
reference genome, by changing **inputs and config only, with no code edits**.

Everything upstream of physical oligo synthesis is in scope. Synthesis is outsourced, so the final
deliverables are **order-ready sequence files** (`baits.fasta` / `baits.csv`), a QC dashboard, and a
clean **`final_covered_targets.bed`** handoff for a downstream variant-analysis pipeline (separate task).

> ⚠️ CaptureForge **designs** baits. It does not order, synthesize, or contact any vendor, and it does
> not run downstream variant analysis.

---

## Pipeline at a glance

```
gene list ─▶ [1] target definition ─▶ [2] masking ─▶ [3] candidate tiling ─▶
[4] specificity ─▶ [5] composition ─▶ [6] gap re-tiling ─▶ [7] boost/dedup ─▶
[8] synthesis scaffolding ─▶ [9] QC aggregation ─▶ [10] order package
```

Each stage is a Nextflow (DSL2) module with QC outputs and explicit BED/FASTA artifacts. Every target
interval carries a **design class** (`coding`, `promoter`, `regulatory`, `cnv`, `str`, `anchor`) that
drives tiling density, masking policy, GC policy, and boosting.

---

## Install

```bash
# 1. Nextflow (needs Java 17+)
curl -s https://get.nextflow.io | bash && sudo mv nextflow /usr/local/bin/

# 2. conda/mamba environments (pinned) — or use the container (see env/)
mamba env create -f env/captureforge.yml      # creates the `captureforge` env

# 3. (optional) containers
docker build -t captureforge:0.1.0 -f env/Dockerfile .        # Mac dev
apptainer build captureforge.sif env/captureforge.def         # HPC
```

## Quickstart

### ▶ Starting from **just a gene list** (recommended on-ramp for new panels)

You have nothing but gene symbols, one per line, in `my_genes.txt`:

```bash
# 1. Scaffold an editable design-spec CSV pre-filled with default design classes
python bin/init_gene_list.py my_genes.txt -o my_panel.spec.csv

# 2. (optional) open my_panel.spec.csv and add promoter / cnv / str / anchor intent per gene

# 3. Build the specificity index FROM THE NO-ALT genome (one-time; required)
bin/build_specificity_index.sh --fasta refs/GRCh38.noalt.fa --aligner blast --out refs/grch38_noalt

# 4. Run the pipeline
nextflow run main.nf -profile mac_local,conda \
    --input my_panel.spec.csv \
    --genome_fasta refs/GRCh38.noalt.fa --gtf refs/GRCh38.112.gtf.gz \
    --specificity_db $PWD/refs/grch38_noalt \
    --outdir results
```

A bare `my_genes.txt` works directly too (applies the default `coding + promoter` classes):
`--input my_genes.txt`.

> **Specificity reference must be the NO-ALT primary assembly == the design genome.** The pipeline
> *enforces* this (fails loudly on alt/decoy/HLA contigs) — a with-alt/toplevel index makes MHC-region
> genes (e.g. GENE) look non-specific against their own alt haplotypes. `--specificity_db` must be an
> **absolute** path to the index prefix (it is referenced in place, not staged).

### ▶ Test profile (10 genes, runs on a Mac in ~2 min)

```bash
bash test/data/make_test_reference.sh        # build a tiny mini-genome fixture from your local GRCh38
nextflow run main.nf -profile test,conda     # or: -profile test (see Apple-Silicon note below)
bash test/check_outputs.sh                    # sanity-checks; expects all 6 design classes
```

**Apple Silicon (arm64) note:** the pinned env builds under Rosetta (osx-64). With `-profile test,conda`
Nextflow manages the env. If instead you pre-activate the env and run `-profile test` (no conda engine),
ensure the env's bin is first so `/usr/bin/env python3` resolves to it:
`conda activate captureforge && export PATH="$CONDA_PREFIX/bin:$PATH"`. The container path (`-profile
test,docker`) avoids this entirely.

### ▶ Full panel on the HPC (SLURM + Apptainer)

See **`docs/hpc_deployment.md`** for the full central-install recipe. Short form:

```bash
module load captureforge/0.1.0           # or: nextflow run /shared/apps/captureforge/main.nf
captureforge-run -params-file my_panel.params.yaml \
    --genome_fasta   $CF_REFS/GRCh38.primary_assembly.fa \
    --gtf            $CF_REFS/GRCh38.112.gtf.gz \
    --specificity_db $CF_REFS/grch38_noalt_bt2 --specificity_aligner bowtie2 \
    --slurm_partition <queue>
```

---

## Input on-ramps

CaptureForge auto-detects and validates any **one** of these (schema: `assets/gene_list_schema.md`):

| Input | Format | Use |
|---|---|---|
| **Plain text** | one gene symbol per line | zero-friction; applies default classes |
| **CSV/TSV spec** | header `gene` + optional `promoter,regulatory,cnv,str,extra_regions,anchors` | clean reusable format |
| **Excel workbook** | curated columns (`Gene`, `Variant class`, `Key variants to capture`, `Locus`, `Burden group`) | rich curated panels |

The logic that turns genes into labelled target regions is **Stage 1 of the pipeline itself** — no prior
project file or script is ever required.

---

## Reference data

All reference paths are params (default GRCh38, overridable for any build/organism). Heavy downloads are
**never run automatically**: `bin/fetch_references.sh` fetches and checksums each reference and writes a
provenance file. See `docs/design_decisions.md` for exactly which releases are used.

---

## Outputs (`results/`)

| Artifact | What |
|---|---|
| `baits.fasta`, `baits.csv` | order-ready bait sequences + metadata |
| `final_covered_targets.bed` | **handoff** target file for the downstream analysis pipeline |
| `qc_dashboard.html` | self-contained visual QC report |
| `figures/*.png,*.svg` | per-gene coverage tracks, GC/Tm/off-target/copy-number plots |
| `metrics.json`, `*.tsv` | machine-readable QC metrics |
| `MANIFEST`, `provenance.json` | md5 checksums + tool/reference/param provenance |

---

## Project layout

See `docs/design_decisions.md` for the running log of every default and its rationale, and the
in-repo directory tree described in the build brief.

## License

MIT (placeholder — see `LICENSE`).
