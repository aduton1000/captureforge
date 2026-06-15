"""captureforge_cli — the unified `captureforge` command-line front door.

A thin wrapper around the existing CaptureForge entry points (the Nextflow `main.nf`
design pipeline, `bin/init_gene_list.py`, `bin/build_specificity_index.sh`). It does NOT
re-architect the pipeline and does NOT change stage logic — it just launches the real
commands so users type `captureforge run/panel/index …` instead of the bare nextflow
form + the conda-PATH dance. The raw forms keep working unchanged.
"""
__version__ = "0.1.0"
