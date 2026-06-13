#!/usr/bin/env bash
# Lightweight sanity check for the `test` profile outputs.
#   nextflow run main.nf -profile test,conda && bash test/check_outputs.sh
set -uo pipefail
OUT="${1:-results}"
fail=0
chk() { if eval "$2"; then echo "  ✓ $1"; else echo "  ✗ $1"; fail=1; fi; }

echo "CaptureForge test-profile sanity check (outdir=$OUT)"
chk "order package: baits.fasta exists & non-empty" "[ -s $OUT/order_package/baits.fasta ]"
chk "order package: baits.csv exists"               "[ -s $OUT/order_package/baits.csv ]"
chk "handoff: final_covered_targets.bed exists"     "[ -s $OUT/stage9_qc/final_covered_targets.bed ]"
chk "metrics.json exists"                           "[ -s $OUT/stage9_qc/metrics.json ]"
chk "MANIFEST exists"                               "[ -s $OUT/order_package/MANIFEST ]"
chk "QC dashboard exists"                           "[ -s $OUT/qc_dashboard.html ]"
chk "figures present (≥8 PNG)"                       "[ \$(ls $OUT/figures/*.png 2>/dev/null | wc -l) -ge 8 ]"

# all six design classes should appear among the candidate/target tags
classes=$(cut -f4 $OUT/stage1_targets/targets.raw.bed 2>/dev/null | awk -F'|' '{print $2}' | sort -u | tr '\n' ' ')
echo "  design classes seen in targets: ${classes:-<none>}"
for c in coding promoter regulatory cnv str anchor; do
  chk "class present: $c" "echo '$classes' | grep -qw $c"
done

# baits should be present in the order package
nb=$(grep -c '^>' "$OUT/order_package/baits.fasta" 2>/dev/null); nb=${nb:-0}
chk "≥1 final bait exported (got $nb)" "[ \"$nb\" -ge 1 ]"

# baits.csv has the expected header columns
hdr=$(head -1 $OUT/order_package/baits.csv 2>/dev/null)
chk "baits.csv header has off_target_flag" "echo '$hdr' | grep -q off_target_flag"

echo
if [ $fail -eq 0 ]; then echo "ALL CHECKS PASSED ✓"; else echo "SOME CHECKS FAILED ✗"; fi
exit $fail
