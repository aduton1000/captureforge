# CaptureForge — HPC deployment (SLURM + Apptainer/Singularity)

A large whole-genome production run is best suited to the HPC: whole-genome indexing and specificity mapping
are heavy. This is a single central install other users invoke.

## 0. One-time: central install

```bash
# clone into a shared location (read-execute for the group)
cd /shared/apps
git clone <your-private-remote>/captureforge.git
cd captureforge

# pin the pipeline version
git checkout v0.1.0     # or a tagged release
```

## 1. Build the container image (native linux-64 toolchain; Docker usually banned on HPC)

The same `env/captureforge.def` builds one `.sif` that runs under **either** Apptainer **or**
SingularityCE — use whichever engine your cluster provides:

```bash
cd /shared/apps/captureforge

# Apptainer sites:
apptainer build captureforge.sif env/captureforge.def

# SingularityCE sites (no `apptainer` binary):
sudo singularity build captureforge.sif env/captureforge.def
# unprivileged alternative if you lack sudo / fakeroot is unavailable:
#   singularity build --fakeroot captureforge.sif env/captureforge.def

# → contains bedtools, samtools, bowtie2, bwa-mem2, genmap, blast(dustmasker),
#   primer3, pysam, Biopython, pyfaidx, pandas, jinja2, matplotlib, viennarna …
```

## 2. Centralized reference data (fetch once, shared read-only)

```bash
export CF_REFS=/shared/refs/captureforge/GRCh38
mkdir -p "$CF_REFS"
# Genome (NO-ALT primary assembly — required) + GTF + masks + dbSNP
bin/fetch_references.sh --outdir "$CF_REFS" genome gtf rmsk simplerepeat segdups dbsnp
# Decompress + index the genome, then build the specificity index FROM THE NO-ALT genome:
gunzip -k "$CF_REFS"/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx "$CF_REFS"/Homo_sapiens.GRCh38.dna.primary_assembly.fa
bin/build_specificity_index.sh \
    --fasta "$CF_REFS"/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --aligner bowtie2 --out "$CF_REFS"/grch38_noalt_bt2
# (blast db alternative: --aligner blast --out "$CF_REFS"/grch38_noalt)
```

> ⚠ The pipeline **enforces** that the specificity reference is no-alt and matches the design genome
> (`specificity_filter.py::assert_noalt_reference`). Pointing it at a toplevel/with-alt index fails loudly.

## 3. Resolve rsID anchors (lightweight; or use the dbSNP VCF)

```bash
bin/resolve_rsids.py --from-targets <(bin/parse_gene_list.py --input <panel> --gtf "$CF_REFS"/*.gtf.gz -o /dev/stdout) \
    -o "$CF_REFS"/anchors.rsid.bed --manifest "$CF_REFS"/refs.provenance.jsonl
# OR pass --dbsnp_vcf "$CF_REFS"/homo_sapiens.dbsnp.GRCh38.vcf.gz to the pipeline instead.
```

## 4. Run on SLURM

Use the same `captureforge run` front door as a local run, pointing `--pipeline` at the installed `main.nf`
(the explicit form, before the Lmod module is loaded — §5 covers the module-based shortcut). Compose the
packaging profile that matches your engine: **`hpc_slurm,singularity`** on a SingularityCE cluster, or
**`hpc_slurm,apptainer`** on an Apptainer one (both run the same `.sif`):

```bash
cd /scratch/$USER/run1
captureforge run \
    --pipeline /shared/apps/captureforge/main.nf \
    -profile hpc_slurm,singularity \
    -params-file my_panel.params.yaml \
    --genome_fasta   $CF_REFS/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --gtf            $CF_REFS/Homo_sapiens.GRCh38.112.gtf.gz \
    --specificity_db $CF_REFS/grch38_noalt_bt2 --specificity_aligner bowtie2 \
    --rmsk_bed       $CF_REFS/rmsk.GRCh38.txt.gz \
    --simplerepeat_bed $CF_REFS/simpleRepeat.GRCh38.txt.gz \
    --rsid_bed       $CF_REFS/anchors.rsid.bed \
    --slurm_partition <queue> --scratch_dir $TMPDIR \
    --outdir results
```

Tune in `conf/hpc_slurm.config`: `slurm_partition` (your queue — e.g. `global` on some clusters, `normal` is
the placeholder default; always overridable with `--slurm_partition <queue>`), `slurm_queue_opts`,
`scratch_dir`, `max_cpus/memory/time`, `container_cache` (the shared image cache; `apptainer_cache` remains a
back-compat alias). Heavy steps (`label index` = specificity, re-tiling) get a large slot (16 cpu / 64 GB).

## 5. Shared invocation: wrapper on PATH (no module system) or Lmod

`share/bin/captureforge-run` is a thin launcher that forwards to `nextflow run … -profile hpc_slurm,<engine>`.
It defaults to the **singularity** engine and resolves the install root, `.sif`, and refs from environment
variables (with install-relative defaults) — so it needs **no module system**.

**No-module clusters (no Lmod):** put the wrapper's directory on a shared `PATH` and run it directly:

```bash
export PATH=/shared/apps/captureforge/share/bin:$PATH
# optional overrides (otherwise derived from the install location):
export CAPTUREFORGE_SIF=/shared/apps/captureforge/captureforge.sif
export CAPTUREFORGE_REFS=/shared/refs/captureforge/GRCh38
# engine defaults to singularity; an Apptainer site sets:
#   export CAPTUREFORGE_ENGINE=apptainer
captureforge-run -params-file my_panel.params.yaml --slurm_partition <queue>
```

**Lmod sites:** install the modulefile stub so users `module load captureforge` (it sets the same env vars),
then run the wrapper:

```bash
module use /shared/apps/captureforge/share/modulefiles
module load captureforge/0.1.0
captureforge-run -params-file my_panel.params.yaml --slurm_partition <queue>
```
