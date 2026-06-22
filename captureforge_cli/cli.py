#!/usr/bin/env python3
"""captureforge — unified command-line front door for the CaptureForge design pipeline.

Subcommands:
  captureforge run     run the design pipeline (thin wrapper over `nextflow run main.nf`)
  captureforge panel   scaffold an editable panel-spec CSV (bin/init_gene_list.py)
  captureforge index   build the Stage-4 specificity index (bin/build_specificity_index.sh)
  captureforge --version | --help

Design notes
------------
* FRIENDLY FRONT DOOR, not a re-architecture. The raw forms keep working unchanged:
  `./.tools/nextflow run main.nf …`, `nextflow run main.nf …`, `python bin/init_gene_list.py …`,
  `bin/build_specificity_index.sh …`.
* CaptureForge is a single-purpose design pipeline (one unnamed `workflow {}`, no `-entry`
  points) — so there are NO per-stage subcommands, by design.
* Nextflow resolution: PREFER the vendored, version-pinned `./.tools/nextflow`; fall back
  to a system `nextflow` on PATH. The resolved launcher is printed for transparency.
"""
import argparse
import importlib.util
import os
import shlex
import shutil
import subprocess
import sys

from . import __version__

PKG_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(PKG_DIR)          # captureforge_cli/ lives at the repo root
BIN = os.path.join(REPO_ROOT, "bin")
VENDORED_NF = os.path.join(REPO_ROOT, ".tools", "nextflow")


# ----------------------------------------------------------------- helpers
def _resolve_nextflow():
    """Return (path, kind). PREFER the vendored launcher (pinned/tested), then system."""
    if os.path.isfile(VENDORED_NF) and os.access(VENDORED_NF, os.X_OK):
        return VENDORED_NF, "vendored (.tools/nextflow)"
    sysnf = shutil.which("nextflow")
    if sysnf:
        return sysnf, "system PATH"
    return None, "none"


def _load_bin(script):
    path = os.path.join(BIN, script)
    if not os.path.isfile(path):
        sys.stderr.write(f"[captureforge] ERROR: cannot find {path}\n")
        sys.exit(2)
    spec = importlib.util.spec_from_file_location(script.replace(".py", ""), path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _run_state(passthrough):
    """Pin Nextflow's run state to a USER-WRITABLE dir, never the read-only shared
    install. Nextflow anchors .nextflow/ (history, history.lock, cache) and the log to
    launchDir == the process CWD, and -work-dir moves ONLY work/ — not those. So we
    launch from the user's CWD (never cd into the install) and, belt-and-suspenders,
    point the log + work dir + the launcher's own scratch (NXF_HOME) at user space.

    Returns (rundir, core_flags, run_flags, env): core_flags go BEFORE `run` (-log is a
    core option), run_flags after it (-work-dir is a run option)."""
    rundir = os.environ.get("CAPTUREFORGE_RUNDIR") or os.getcwd()
    core, run = [], []
    has_log = any(a == "-log" or a.startswith("-log=") for a in passthrough)
    has_work = any(a in ("-work-dir", "-w") or a.startswith("-work-dir=") for a in passthrough)
    if not has_log:
        core += ["-log", os.path.join(rundir, ".nextflow.log")]
    if not has_work and not os.environ.get("NXF_WORK"):
        run += ["-work-dir", os.path.join(rundir, "work")]
    env = os.environ.copy()
    env.setdefault("NXF_HOME", os.path.join(os.path.expanduser("~"), ".nextflow"))
    return rundir, core, run, env


def _require_writable(rundir):
    """Refuse to run from a non-writable dir with a clear message — instead of letting
    Nextflow emit the cryptic '.nextflow/history.lock (Permission denied)'."""
    if not os.access(rundir, os.W_OK):
        sys.stderr.write(
            f"[captureforge] ERROR: run directory is not writable: {rundir}\n"
            "  CaptureForge writes run state (.nextflow/, logs, work/) to this directory.\n"
            "  Run from a writable directory (your home or a project dir), e.g.:\n"
            "    cd ~ && captureforge run -params-file site.params.yaml --input panel.txt --outdir results\n"
            "  (or set CAPTUREFORGE_RUNDIR=/path/to/writable to override).\n")
        sys.exit(1)


def _warn_env_if_local(passthrough):
    """arm64 backstop: warn (do NOT auto-activate) if the captureforge conda env / its
    tools are not on PATH for a LOCAL (non-container) run. The user stays in control.

    Suppressed for ANY container engine (the tools come from the image, not a host
    conda env) — including singularity, which the old check missed and so emitted a
    misleading 'conda activate captureforge' on every in-container cluster run."""
    joined = " ".join(passthrough)
    # A container engine selected via -profile ...,<engine>, an explicit -with-<engine>,
    # or a wired --container_image all mean the tools live in the image — no conda needed.
    engines = ("docker", "apptainer", "singularity", "podman", "shifter", "charliecloud", "sarus")
    if any(e in joined for e in engines) or "--container_image" in passthrough:
        return                                  # a container provides the tools
    env = os.environ.get("CONDA_DEFAULT_ENV", "")
    prefix = os.environ.get("CONDA_PREFIX", "")
    if env == "captureforge" or prefix.rstrip("/").endswith("captureforge"):
        return
    sys.stderr.write(
        "[captureforge] WARN: the `captureforge` conda env was not detected. Local "
        "(non-container) runs — especially Apple-Silicon arm64 — need its tools on PATH:\n"
        "    conda activate captureforge\n"
        "  Proceeding anyway (use a container engine, e.g. -profile <exec>,docker, to avoid this).\n")


# Pipeline parameters are declared in Nextflow (nextflow.config / conf/), NOT in
# argparse — `captureforge run` forwards them verbatim. argparse therefore can't own
# their --help, so we document them here (printed as the `run` epilog) WITH provenance
# so a reader knows, per flag, whether they set it or it is already wired.
#   [install]         once in the site env file (CAPTUREFORGE_HOME/_SIF/_REFS/_BIND/_ENGINE)
#   [site]            once per cluster in the params file (-params-file)
#   [per-run]         you provide every run (the minimal set)
#   [optional]        has a sane default; set only to override (default shown)
#   [class-dependent] only needed if the panel uses that design class
# Defaults below are the ACTUAL nextflow.config / conf defaults (no invented values).
PASSTHROUGH_HELP = """\
Nextflow pipeline parameters (forwarded to `nextflow run`; declared in
nextflow.config / conf/, so they are passthrough — not listed under "options" above).
Provenance: [install]=env file once · [site]=params file once per cluster ·
[per-run]=every run · [optional]=defaulted override · [class-dependent]=only for that class.

  -profile NAMES        [per-run] execution+packaging, composed: test|local|hpc_slurm
                        with conda|docker|apptainer|singularity (e.g. local,conda)
  -params-file FILE     [site] YAML carrying the site values (genome, index, binds, …)

 references
  --input PATH          [per-run] REQUIRED gene list / CSV-TSV spec / .xlsx (auto-detected)
  --genome_fasta PATH   [site] REQUIRED no-alt primary-assembly FASTA (cleaned headers)
  --gtf PATH            [site] REQUIRED Ensembl GTF (gene models)
  --specificity_db PFX  [site] aligner index / BLAST db prefix (absolute, no-alt). No
                        hard default — Stage 4/6 fail without it (main.nf warns)
  --specificity_aligner X [site] bowtie2 (default) | bwa | blast
  --genome_build STR    [optional] default GRCh38 (label only)
  --ensembl_release N   [optional] default 112
  --mane_gff PATH       [optional] MANE Select GFF3 for canonical pick (default: none)
  --regulatory_gff PATH [class-dependent] Ensembl Regulatory Build — regulatory class
  --rmsk_bed PATH       [class-dependent] UCSC RepeatMasker — Stage-2 masking / str class
  --simplerepeat_bed PATH [class-dependent] UCSC simpleRepeat — STR-core removal
  --segdups_bed PATH    [class-dependent] UCSC genomicSuperDups — cnv masking
  --dbsnp_vcf PATH      [class-dependent] dbSNP VCF for rsID anchors (anchor class)
  --rsid_bed PATH       [class-dependent] pre-resolved rsID BED (anchor class)
                        (optional refs absent → that feature is skipped, no error)

 container / HPC
  --container_image X   [site] image every task runs in (.sif path | docker tag).
                        REQUIRED when a container engine is enabled (else fails loudly);
                        the wrapper supplies it from CAPTUREFORGE_SIF
  --bind_paths PATHS    [site] host path(s) to bind into containers for in-place refs
                        (e.g. the index). Default: none. Comma-separate multiple roots
  --slurm_partition Q   [site] SLURM queue. REQUIRED under -profile hpc_slurm (no safe
                        default; fails loudly if unset). Find one with `sinfo -s`
  --slurm_queue_opts S  [optional] extra sbatch/clusterOptions (e.g. --account=...)
  --scratch_dir PATH    [optional] per-task scratch under SLURM (default /tmp)

 design thresholds (all [optional]; defaults shown)
  --bait_len N (120)  --splice_flank N (15)  --canonical_source X (mane)
  --promoter_upstream N (2000)  --promoter_downstream N (200)  --anchor_flank N (150)
  --gc_min F (0.25)  --gc_max F (0.70)  --gc_max_promoter F (0.80)
  --tm_min F (60.0)  --tm_max F (90.0)  --homopolymer_max N (10)
  --offtarget_min_identity F (0.90)  --offtarget_min_lenfrac F (0.50)  --offtarget_max_hits N (1)
  --retile_max_passes N (3)  --coverage_target F (0.95)
  --bait_format X (ivt_rna: dna_pool|ivt_rna|ready_to_use)  --synth_max_len N (170)
  --per_oligo_price F (0.00200)  --pool_tiers LIST ([12000,92000,244000])

 output
  --panel_name STR      [per-run/optional] label on outputs (default captureforge_panel)
  --outdir PATH         [per-run] results directory (default results)

See `params.example.yaml`, `conf/cluster.params.example.yaml`, and the manual's
Parameter Reference for the full table. Prefers the vendored .tools/nextflow;
falls back to a system nextflow on PATH."""


# ----------------------------------------------------------------- subcommands
def cmd_run(rest):
    ap = argparse.ArgumentParser(
        prog="captureforge run", allow_abbrev=False, add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Run the CaptureForge design pipeline. Thin wrapper over "
            "`nextflow run main.nf`; -profile, -params-file, -resume and every "
            "--<param> pass straight through to Nextflow.\n\n"
            "Most options are set ONCE (env file + params file), so a typical run is short:\n"
            "    captureforge run -profile local,conda -params-file params.yaml \\\n"
            "        --input panel.txt --outdir results\n"
            "On a deployed cluster the wrapper makes it shorter still:\n"
            "    captureforge-run -params-file site.params.yaml --input panel.txt --outdir results\n"
            "The flags below OVERRIDE params-file values; you rarely type them all."),
        epilog=PASSTHROUGH_HELP)
    ap.add_argument("--pipeline", metavar="PATH",
                    help="[optional] pipeline script to run "
                         "(default: this repo's main.nf; point at the installed "
                         "main.nf for an out-of-tree deploy)")
    ap.add_argument("--print-cmd", action="store_true",
                    help="[optional] print the resolved nextflow command and exit "
                         "(do not execute) — same as setting CAPTUREFORGE_DRY_RUN")
    known, passthrough = ap.parse_known_args(rest)

    nf, kind = _resolve_nextflow()
    if nf is None:
        sys.stderr.write("[captureforge] ERROR: no Nextflow found — the vendored "
                         ".tools/nextflow is missing/not executable and no system "
                         "`nextflow` is on PATH. Install Nextflow (Java 17+).\n")
        return 127
    pipeline = known.pipeline or os.path.join(REPO_ROOT, "main.nf")
    if not known.pipeline and not os.path.isfile(pipeline):
        pipeline = "main.nf"
    # Run from the user's CWD (launchDir), NOT the install — main.nf is referenced by
    # absolute path so .nextflow/ + logs + work/ land in user space, not the repo.
    rundir, core, run_flags, env = _run_state(passthrough)
    cmd = [nf, *core, "run", pipeline, *run_flags, *passthrough]

    sys.stderr.write(f"[captureforge] nextflow resolved: {kind} -> {nf}\n")
    sys.stderr.write("[captureforge] + " + " ".join(shlex.quote(c) for c in cmd) + "\n")
    if known.print_cmd or os.environ.get("CAPTUREFORGE_DRY_RUN"):
        print(" ".join(shlex.quote(c) for c in cmd))
        return 0
    _warn_env_if_local(passthrough)
    _require_writable(rundir)
    try:
        return subprocess.call(cmd, env=env)
    except FileNotFoundError:
        sys.stderr.write(f"[captureforge] ERROR: could not execute {nf}\n")
        return 127


def cmd_panel(rest):
    # scaffold the panel-spec CSV. init_gene_list.main(argv) is importable; the raw form
    # `python bin/init_gene_list.py …` still works.
    mod = _load_bin("init_gene_list.py")
    return mod.main(rest) or 0


def cmd_index(rest):
    # build the Stage-4 specificity index (thin passthrough to the real bash script).
    script = os.path.join(BIN, "build_specificity_index.sh")
    if not os.path.isfile(script):
        sys.stderr.write(f"[captureforge] ERROR: cannot find {script}\n")
        return 2
    cmd = [script, *rest]
    sys.stderr.write("[captureforge] + " + " ".join(shlex.quote(c) for c in cmd) + "\n")
    try:
        return subprocess.call(cmd, cwd=REPO_ROOT)
    except FileNotFoundError:
        sys.stderr.write(f"[captureforge] ERROR: could not execute {script}\n")
        return 127


# ----------------------------------------------------------------- dispatch
SUBCOMMANDS = {
    "run": (cmd_run, "run the design pipeline (wrapper over `nextflow run main.nf`)"),
    "panel": (cmd_panel, "scaffold an editable panel-spec CSV from gene symbols"),
    "index": (cmd_index, "build the Stage-4 specificity index from the design genome"),
}


def _top_help():
    lines = ["captureforge — front door for the CaptureForge bait-design pipeline",
             f"version {__version__}", "",
             "usage: captureforge <command> [options]", "", "commands:"]
    for name, (_, desc) in SUBCOMMANDS.items():
        lines.append(f"  {name:<7} {desc}")
    lines += ["", "  --version    print version", "  --help       show this help", "",
              "Raw forms still work: `./.tools/nextflow run main.nf …`, "
              "`python bin/init_gene_list.py …`, `bin/build_specificity_index.sh …`.",
              "Run `captureforge <command> --help` for command-specific options."]
    return "\n".join(lines)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help"):
        print(_top_help())
        return 0
    if argv[0] in ("-V", "--version"):
        print(__version__)
        return 0
    cmd, rest = argv[0], argv[1:]
    if cmd not in SUBCOMMANDS:
        sys.stderr.write(f"[captureforge] unknown command '{cmd}'\n\n{_top_help()}\n")
        return 2
    return SUBCOMMANDS[cmd][0](rest)


if __name__ == "__main__":
    sys.exit(main())
