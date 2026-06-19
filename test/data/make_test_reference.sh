#!/usr/bin/env bash
# Rebuild the tiny generic mini-genome fixture for the `test` profile.
#
# NOTE: the test profile ships a PRE-BUILT fixture (test/data/ref/*) — the
# quickstart `captureforge run -profile test` does NOT need this script. It is
# provided only for reproducibility / to refresh the fixture after editing
# test_genes.csv.
#
# It carves each gene in test_genes.csv (± flank) out of a full GRCh38 (no-alt)
# FASTA + Ensembl GTF into its own local contig, so the full 10-stage pipeline
# runs end-to-end in minutes. The genes are generic, non-panel exemplars
# (TP53/EGFR/BRCA1/…) and the anchors are synthetic placeholders shipped in
# test/data/anchors.generic.bed — no project panel data is referenced.
#
# Point it at your own local references (it does NOT assume any machine path):
#   FULL_FASTA=/path/GRCh38.fa FULL_GTF=/path/Homo_sapiens.GRCh38.gtf.gz \
#       bash test/data/make_test_reference.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"

# Required: point these at YOUR local GRCh38 no-alt FASTA + Ensembl GTF.
FULL_FASTA="${FULL_FASTA:-}"
FULL_GTF="${FULL_GTF:-$ROOT/refs/Homo_sapiens.GRCh38.112.gtf.gz}"
# Default to the shipped generic anchor placeholders (not project panel anchors).
RSID_BED="${RSID_BED:-$HERE/anchors.generic.bed}"
FLANK="${FLANK:-30000}"

[[ -n "$FULL_FASTA" && -s "$FULL_FASTA" ]] || {
    echo "ERROR: set FULL_FASTA=/path/to/GRCh38(no-alt).fa  (Ensembl-style contig names: 1,2,…,X)" >&2; exit 1; }
[[ -s "$FULL_GTF" ]] || { echo "ERROR: set FULL_GTF=/path/to/Ensembl.gtf.gz" >&2; exit 1; }

rsid_arg=()
[[ -s "$RSID_BED" ]] && rsid_arg=(--rsid-bed "$RSID_BED")

python "$ROOT/bin/make_test_reference.py" \
    --genes "$HERE/test_genes.csv" \
    --gtf "$FULL_GTF" --fasta "$FULL_FASTA" \
    "${rsid_arg[@]}" --flank "$FLANK" --outdir "$HERE/ref"
echo "[make_test_reference] done → $HERE/ref"
