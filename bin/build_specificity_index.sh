#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════════
#  CaptureForge — bin/build_specificity_index.sh
#  Build the Stage-4 specificity index FROM THE DESIGN GENOME so specificity is
#  measured against exactly the reference the baits target.
#
#  IMPORTANT: use the NO-ALT primary assembly. Indexing a toplevel/with-alt
#  reference makes genes in the MHC and other alt-haplotype regions look
#  non-specific against their own alternate-haplotype copies (false off-targets).
#
#  Provenance: a one-off [site] setup step, run ONCE per cluster/genome (not per run).
#  The prefix it writes becomes the `specificity_db` param in your params file.
#
#  Usage:
#    bin/build_specificity_index.sh --fasta GRCh38.primary.fa --aligner blast   --out refs/grch38_noalt
#    bin/build_specificity_index.sh --fasta GRCh38.primary.fa --aligner bowtie2 --out refs/grch38_bt2
#  Options:
#    --fasta PATH    [required] no-alt primary-assembly FASTA (same cleaned FASTA as genome_fasta)
#    --aligner NAME  [optional] blast (default) | bowtie2 | bwa
#    --out PREFIX    [required] output index/db prefix (→ set as specificity_db)
#    THREADS=N (env) [optional] threads for bowtie2-build (default 4)
# ═══════════════════════════════════════════════════════════════════════════
set -euo pipefail
FASTA=""; ALIGNER="blast"; OUT=""
while [[ $# -gt 0 ]]; do case "$1" in
  --fasta) FASTA="$2"; shift 2;;
  --aligner) ALIGNER="$2"; shift 2;;
  --out) OUT="$2"; shift 2;;
  -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
  *) echo "unknown arg: $1" >&2; exit 1;;
esac; done
[[ -s "$FASTA" ]] || { echo "ERROR: --fasta required and must exist" >&2; exit 1; }
[[ -n "$OUT" ]] || { echo "ERROR: --out prefix required" >&2; exit 1; }
mkdir -p "$(dirname "$OUT")"

# warn if the FASTA looks like it contains alt contigs
if grep -qiE 'HSCHR|_alt|MHC_' "${FASTA}.fai" 2>/dev/null || grep -qiE '^>.*(HSCHR|_alt|MHC_)' "$FASTA"; then
  echo "WARNING: '$FASTA' appears to contain ALT contigs — specificity will over-reject MHC/segdup genes." >&2
fi

# Warn if FASTA headers are NOT first-token-clean. Index builders store the FULL
# header as each contig name, but `samtools faidx` keeps only the first token —
# a mismatch makes Stage 4's same-reference check see two different references.
# (Stage 4 now also compares on first-tokens defensively, but cleaning is best.)
if grep -m1 -qE '^>[^[:space:]]+[[:space:]]+\S' "$FASTA"; then
  echo "WARNING: '$FASTA' has multi-token FASTA headers; the index will store the full" >&2
  echo "         header as each contig name while samtools faidx keeps only the first token." >&2
  echo "         Clean headers first (and re-faidx) so genome and index names match:" >&2
  echo "             awk '/^>/{print \$1; next} {print}' in.fa > clean.fa && samtools faidx clean.fa" >&2
  echo "         (bin/fetch_references.sh genome writes a *.clean.fa automatically.)" >&2
fi

case "$ALIGNER" in
  blast)    echo "[build] makeblastdb → ${OUT}";   makeblastdb -in "$FASTA" -dbtype nucl -out "$OUT" -title "$(basename "$OUT")" ;;
  bowtie2)  echo "[build] bowtie2-build → ${OUT}";  bowtie2-build --threads "${THREADS:-4}" "$FASTA" "$OUT" ;;
  bwa)      echo "[build] bwa index → ${OUT}";      bwa index -p "$OUT" "$FASTA" ;;
  *) echo "ERROR: --aligner must be blast|bowtie2|bwa" >&2; exit 1;;
esac
echo "[build] done: ${OUT} (aligner=${ALIGNER})"
