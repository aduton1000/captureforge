#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════════
#  CaptureForge — bin/fetch_references.sh
#  Download + checksum reference data, write a provenance file. SELECTIVE and
#  NOT auto-run by the pipeline: heavy downloads happen only when you ask.
#
#  Usage:
#     bin/fetch_references.sh --outdir refs --release 112 [targets...]
#  Targets (default: gtf):
#     gtf regulatory rmsk segdups dbsnp genome  | all
#  Examples:
#     bin/fetch_references.sh --outdir refs gtf            # just the GTF (~50 MB)
#     bin/fetch_references.sh --outdir refs all            # everything (heavy!)
#     bin/fetch_references.sh --record-local /path/ref.fa  # md5 an existing FASTA
# ═══════════════════════════════════════════════════════════════════════════
set -euo pipefail

OUTDIR="refs"
RELEASE="112"
BUILD="GRCh38"
SPECIES="homo_sapiens"
RECORD_LOCAL=""
TARGETS=()

err()  { echo "ERROR: $*" >&2; exit 1; }
log()  { echo "[fetch_references] $*" >&2; }

# md5 wrapper (macOS `md5 -q` vs linux `md5sum`)
md5of() { if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | awk '{print $1}'; else md5 -q "$1"; fi; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir)       OUTDIR="$2"; shift 2;;
    --release)      RELEASE="$2"; shift 2;;
    --build)        BUILD="$2"; shift 2;;
    --species)      SPECIES="$2"; shift 2;;
    --record-local) RECORD_LOCAL="$2"; shift 2;;
    -h|--help)      grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
    gtf|regulatory|rmsk|segdups|dbsnp|genome|all) TARGETS+=("$1"); shift;;
    *) err "unknown argument: $1";;
  esac
done
[[ ${#TARGETS[@]} -eq 0 && -z "$RECORD_LOCAL" ]] && TARGETS=("gtf")
[[ " ${TARGETS[*]} " == *" all "* ]] && TARGETS=(genome gtf regulatory rmsk segdups dbsnp)

mkdir -p "$OUTDIR"
PROV="$OUTDIR/refs.provenance.json"
TS="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
ENSEMBL="https://ftp.ensembl.org/pub/release-${RELEASE}"
UCSC="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database"

# provenance accumulates as a JSON array of objects
PROV_ENTRIES=()
record() {  # name url localpath
  local name="$1" url="$2" path="$3" md5; md5="$(md5of "$path")"
  local bytes; bytes="$(stat -f%z "$path" 2>/dev/null || stat -c%s "$path")"
  PROV_ENTRIES+=("{\"name\":\"$name\",\"build\":\"$BUILD\",\"release\":\"$RELEASE\",\"source\":\"$url\",\"path\":\"$path\",\"bytes\":$bytes,\"md5\":\"$md5\",\"fetched\":\"$TS\"}")
  log "✓ $name  md5=$md5  ($bytes bytes)"
}

fetch() {  # url localpath
  local url="$1" path="$2"
  [[ -s "$path" ]] && { log "exists, skip: $path"; return; }
  log "downloading: $url"
  curl -fL --retry 3 -o "$path" "$url" || err "download failed: $url"
}

# ── record an existing local FASTA (no download) ───────────────────────────
if [[ -n "$RECORD_LOCAL" ]]; then
  [[ -s "$RECORD_LOCAL" ]] || err "no such file: $RECORD_LOCAL"
  log "checksumming local reference (may take a moment for a whole genome)…"
  record "genome_fasta_local" "file://$RECORD_LOCAL" "$RECORD_LOCAL"
fi

for t in "${TARGETS[@]:-}"; do case "$t" in
  genome)
    f="$OUTDIR/${SPECIES^}.${BUILD}.dna.primary_assembly.fa.gz"
    fetch "${ENSEMBL}/fasta/${SPECIES}/dna/Homo_sapiens.${BUILD}.dna.primary_assembly.fa.gz" "$f"
    record "genome_fasta" "${ENSEMBL}/fasta/${SPECIES}/dna/" "$f" ;;
  gtf)
    f="$OUTDIR/Homo_sapiens.${BUILD}.${RELEASE}.gtf.gz"
    fetch "${ENSEMBL}/gtf/${SPECIES}/Homo_sapiens.${BUILD}.${RELEASE}.gtf.gz" "$f"
    record "gtf" "${ENSEMBL}/gtf/${SPECIES}/" "$f" ;;
  regulatory)
    f="$OUTDIR/Homo_sapiens.${BUILD}.regulatory_features.gff.gz"
    fetch "${ENSEMBL}/regulation/${SPECIES}/${BUILD}/annotation/Homo_sapiens.${BUILD}.regulatory_features.v${RELEASE}.gff.gz" "$f" \
      || log "regulatory build path varies by release; edit URL if needed"
    [[ -s "$f" ]] && record "regulatory_gff" "${ENSEMBL}/regulation/" "$f" ;;
  rmsk)
    f="$OUTDIR/rmsk.${BUILD}.txt.gz"
    fetch "${UCSC}/rmsk.txt.gz" "$f"
    record "rmsk" "${UCSC}/rmsk.txt.gz" "$f" ;;
  segdups)
    f="$OUTDIR/genomicSuperDups.${BUILD}.txt.gz"
    fetch "${UCSC}/genomicSuperDups.txt.gz" "$f"
    record "segdups" "${UCSC}/genomicSuperDups.txt.gz" "$f" ;;
  dbsnp)
    f="$OUTDIR/dbsnp.${BUILD}.vcf.gz"
    fetch "${ENSEMBL}/variation/vcf/${SPECIES}/${SPECIES}.vcf.gz" "$f"
    record "dbsnp_vcf" "${ENSEMBL}/variation/vcf/${SPECIES}/" "$f" ;;
esac; done

# ── write provenance JSON (merge with any prior entries) ───────────────────
{
  echo "{"
  echo "  \"generated\": \"$TS\","
  echo "  \"build\": \"$BUILD\", \"ensembl_release\": \"$RELEASE\","
  echo "  \"references\": ["
  IFS=$'\n'; first=1
  for e in "${PROV_ENTRIES[@]:-}"; do
    [[ -z "$e" ]] && continue
    [[ $first -eq 0 ]] && echo ","
    printf "    %s" "$e"; first=0
  done
  echo
  echo "  ]"
  echo "}"
} > "$PROV"
log "provenance written → $PROV"
log "DONE. Note these md5s in docs/design_decisions.md and the run report."
