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

> **Shared, read-only conda?** If you install the `captureforge` CLI (or build the env) against a cluster's
> shared conda whose package cache is not group-writable, conda fails with *"cannot set permissions"* on the
> pkgs cache. Point the cache at a writable dir first:
> ```bash
> export CONDA_PKGS_DIRS="$HOME/.conda/pkgs"
> ```
> On a container-only deployment (below) you don't need a host conda env at all — the tools come from the `.sif`.

## 1. Build the container image (native linux-64 toolchain; Docker usually banned on HPC)

The same `env/captureforge.def` builds one `.sif` that runs under **either** Apptainer **or**
SingularityCE — use whichever engine your cluster provides:

```bash
cd /shared/apps/captureforge

# Apptainer sites:
apptainer build captureforge.sif env/captureforge.def

# SingularityCE sites (no `apptainer` binary):
sudo singularity build captureforge.sif env/captureforge.def
# unprivileged alternative (needs fakeroot configured — a /etc/subuid + /etc/subgid
# entry for your user; ask your admin if `--fakeroot` errors):
#   singularity build --fakeroot captureforge.sif env/captureforge.def

# → contains bedtools, samtools, bowtie2, bwa-mem2, genmap, blast(dustmasker),
#   primer3, pysam, Biopython, pyfaidx, pandas, jinja2, matplotlib, viennarna,
#   AND procps (`ps`) + base shell utils — Nextflow calls `ps` to poll tasks, so a
#   stripped image makes every task fail. Build from the provided `.def` (don't trim it).
```

> **Verify a tool is in the image** with a NON-login shell — `singularity exec <img> <tool> --version`. A login
> shell (`bash -lc …`) re-sources profiles and can reset PATH, hiding the image's tools and giving a misleading
> "not found":
> ```bash
> singularity exec captureforge.sif bowtie2 --version     # ✓ what tasks actually see
> singularity exec captureforge.sif python -c 'import matplotlib; print(matplotlib.__version__)'
> singularity exec captureforge.sif ps --version          # procps present (Nextflow needs it)
> ```

## 2. Centralized reference data (fetch once, shared read-only)

```bash
export CF_REFS=/shared/refs/captureforge/GRCh38
mkdir -p "$CF_REFS"
# Genome (NO-ALT primary assembly — required) + GTF + masks + dbSNP.
# The `genome` target also writes a HEADER-CLEANED FASTA (+ .fai):
#   Homo_sapiens.GRCh38.dna.primary_assembly.clean.fa
bin/fetch_references.sh --outdir "$CF_REFS" genome gtf rmsk simplerepeat segdups dbsnp

# Build the specificity index FROM THE CLEANED NO-ALT genome (same file you pass as
# --genome_fasta). bowtie2 is recommended on the HPC:
bin/build_specificity_index.sh \
    --fasta "$CF_REFS"/Homo_sapiens.GRCh38.dna.primary_assembly.clean.fa \
    --aligner bowtie2 --out "$CF_REFS"/grch38_noalt_bt2
# (blast db alternative: --aligner blast --out "$CF_REFS"/grch38_noalt)
```

> **Why the cleaned FASTA (lesson learned).** `bowtie2-build`/`bwa index`/`makeblastdb` store the **full FASTA
> header** as each contig name (e.g. `1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF`), while
> `samtools faidx` keeps only the **first token** (`1`). If the genome and the index disagree on names, Stage 4's
> same-reference check sees two different references and fails. `bin/fetch_references.sh genome` truncates headers
> to the first token (`1 2 … X Y MT` + scaffolds) and faidxes the result — use that `*.clean.fa` for **both**
> `--genome_fasta` and the index build. (Stage 4 also compares names on their first token defensively, so an
> un-cleaned index still works; cleaning at prep time is the supported path. If `samtools` wasn't on PATH at
> fetch time, faidx the cleaned FASTA once inside the container: `singularity exec captureforge.sif samtools
> faidx <clean.fa>`.)

> ⚠ The pipeline **enforces** that the specificity reference is no-alt and matches the design genome
> (`specificity_filter.py::assert_noalt_reference`). Pointing it at a toplevel/with-alt index fails loudly.
> Specificity reporting is **capped** (`bowtie2 -k`, not `-a`) so it scales to the whole genome — a bait in a
> repeat is rejected without enumerating its millions of hits. Nothing to configure; it explains the fast runtime.

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
    --genome_fasta   $CF_REFS/Homo_sapiens.GRCh38.dna.primary_assembly.clean.fa \
    --gtf            $CF_REFS/Homo_sapiens.GRCh38.112.gtf.gz \
    --specificity_db $CF_REFS/grch38_noalt_bt2 --specificity_aligner bowtie2 \
    --rmsk_bed       $CF_REFS/rmsk.GRCh38.txt.gz \
    --simplerepeat_bed $CF_REFS/simpleRepeat.GRCh38.txt.gz \
    --rsid_bed       $CF_REFS/anchors.rsid.bed \
    --slurm_partition <queue> --scratch_dir $TMPDIR \
    --container_image /hpc/opt/captureforge/images/captureforge.sif \
    --outdir results
# refs that live outside $HOME/CWD (e.g. the --specificity_db index under /hpc) must be
# bound into task containers — the reliable way is an exported env var (see callout below):
export SINGULARITY_BIND=/hpc      # APPTAINER_BIND=/hpc on an Apptainer site; [cluster-specific]
```

> **`--container_image` is required when an engine is enabled.** It sets `process.container`, so every task runs
> via `singularity exec <image> …` (the `.sif` you built in §1). Without it, Nextflow turns the engine on but
> wires no image to processes — tasks then run on the **host** (silently, until a stage needs a tool the host
> lacks, e.g. `matplotlib`). The pipeline now **fails loudly** in that case. For singularity/apptainer use a
> local `.sif` path (`/hpc/opt/...` `[cluster-specific]`); for docker use the image tag. Put it in the params
> file as `container_image:` so it's not retyped each run.

> **Container bind path.** Singularity/Apptainer only auto-mount `$HOME`, `/tmp`, and the CWD. The specificity
> index (`--specificity_db`) is **referenced in place — not staged into work dirs** — so if your references live
> outside those (e.g. under `/hpc/opt/...` `[cluster-specific]`) you must bind that root into the task containers,
> or they fail with "No such file or directory" on the index.
>
> - **Primary — `export SINGULARITY_BIND=/hpc`** (or `APPTAINER_BIND` on an Apptainer site). Nextflow scrubs the
>   task environment, but the `hpc_slurm`/`singularity`/`apptainer` profiles **whitelist** these vars
>   (`envWhitelist`) so they survive the scrub and the engine honors them. This is the most robust route.
> - **Alternative — `--bind_paths /hpc`** (or `bind_paths:` in the params file). Applied per task via
>   `process.containerOptions`. (Note: this is *not* wired through engine `runOptions`, which evaluates at
>   config-parse time before the value merges and would never bind — the closure form is what makes it work.)
>
> Either works; you don't need both. Default is empty (no bind). The `captureforge-run` wrapper sets both from
> `CAPTUREFORGE_BIND`.

Tune in `conf/hpc_slurm.config`: `slurm_partition` (**required** for SLURM runs — there is no safe default;
discover yours with `sinfo -s`, e.g. `global` on some clusters, and set `--slurm_partition <queue>` or
`slurm_partition: <queue>` in the params file. If it is unset the run **fails loudly** before any job is
submitted, rather than silently submitting to a non-existent partition), `slurm_queue_opts`,
`scratch_dir`, `max_cpus/memory/time`, `container_cache` (the shared image cache; `apptainer_cache` remains a
back-compat alias). Heavy steps (`label index` = specificity, re-tiling) get a large slot (16 cpu / 64 GB).

## 5. Shared invocation: wrapper on PATH (no module system) or Lmod

`share/bin/captureforge-run` is a thin launcher that forwards to `nextflow run … -profile hpc_slurm,<engine>`.
It defaults to the **singularity** engine and resolves the install root, `.sif`, and refs from environment
variables (with install-relative defaults) — so it needs **no module system**.

**No-module clusters (no Lmod) — source the site env once, then run short.** Copy
`share/captureforge-env.sh.example` to a shared path, fill in the values once, and have users source it from
`~/.bashrc`. After that, `captureforge-run` is on PATH and the bind/image/refs are configured — no per-run
exports:

```bash
# one-time, per site (admin):
cp share/captureforge-env.sh.example /shared/apps/captureforge/captureforge-env.sh
$EDITOR /shared/apps/captureforge/captureforge-env.sh      # set HOME/SIF/REFS/ENGINE/BIND

# one-time, per user:
echo 'source /shared/apps/captureforge/captureforge-env.sh' >> ~/.bashrc

# every run thereafter (queue + refs + image all come from the params file / env):
captureforge-run -params-file site.params.yaml --input panel.txt --outdir results
```

Start `site.params.yaml` from the annotated template **`conf/cluster.params.example.yaml`** (copy and edit
once): it carries `genome_fasta` (the cleaned FASTA), `gtf`, `specificity_db`, `specificity_aligner`,
`bind_paths`, `container_image`, and `slurm_partition` — so they are not retyped per run. The wrapper still
accepts overrides (`--slurm_partition <queue>`, `--outdir …`, any `nextflow` flag).

The env file sets these (all also settable by hand):

```bash
export PATH=/shared/apps/captureforge/share/bin:$PATH   # wrapper on PATH
export CAPTUREFORGE_SIF=/shared/apps/captureforge/captureforge.sif
export CAPTUREFORGE_REFS=/shared/refs/captureforge/GRCh38
export CAPTUREFORGE_ENGINE=singularity                  # apptainer on an Apptainer site
export CAPTUREFORGE_BIND=/hpc                            # [cluster-specific] refs/install root (binds into tasks)
```

**Lmod sites:** install the modulefile stub so users `module load captureforge` (it sets the same env vars),
then run the wrapper:

```bash
module use /shared/apps/captureforge/share/modulefiles
module load captureforge/0.1.0
captureforge-run -params-file my_panel.params.yaml --slurm_partition <queue>
```
