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


# ----------------------------------------------------------------- subcommands
def cmd_run(rest):
    ap = argparse.ArgumentParser(
        prog="captureforge run", allow_abbrev=False, add_help=True,
        description="Run the CaptureForge design pipeline. Thin wrapper over "
                    "`nextflow run main.nf`; all other flags (-profile, -params-file, "
                    "-resume, --<param>) pass straight through to Nextflow.",
        epilog="Prefers the vendored .tools/nextflow; falls back to system nextflow.")
    ap.add_argument("--pipeline", metavar="PATH",
                    help="pipeline script to run (default: the repo's main.nf)")
    ap.add_argument("--print-cmd", action="store_true",
                    help="print the resolved command and exit (do not execute)")
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
    cmd = [nf, "run", pipeline, *passthrough]

    sys.stderr.write(f"[captureforge] nextflow resolved: {kind} -> {nf}\n")
    sys.stderr.write("[captureforge] + " + " ".join(shlex.quote(c) for c in cmd) + "\n")
    if known.print_cmd or os.environ.get("CAPTUREFORGE_DRY_RUN"):
        print(" ".join(shlex.quote(c) for c in cmd))
        return 0
    _warn_env_if_local(passthrough)
    try:
        return subprocess.call(cmd, cwd=REPO_ROOT)
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
