#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════════
#  CaptureForge — bin/fetch_references.sh
#  Download + checksum reference data, MERGING into a provenance manifest. The
#  pipeline never auto-runs this: heavy downloads happen only when you ask.
#
#  Usage:
#     bin/fetch_references.sh --outdir refs --release 112 [targets...]
#  Targets (default: gtf):
#     gtf regulatory rmsk simplerepeat segdups dbsnp genome  | all
#  Examples:
#     bin/fetch_references.sh --outdir refs gtf
#     bin/fetch_references.sh --outdir refs rmsk simplerepeat   # STR/repeat masks
#     bin/fetch_references.sh --outdir refs all                 # everything (heavy!)
#     bin/fetch_references.sh --record-local /path/ref.fa       # md5 an existing FASTA
# ═══════════════════════════════════════════════════════════════════════════
set -euo pipefail

OUTDIR="refs"; RELEASE="112"; BUILD="GRCh38"; SPECIES="homo_sapiens"
RECORD_LOCAL=""; TARGETS=()

err()  { echo "ERROR: $*" >&2; exit 1; }
log()  { echo "[fetch_references] $*" >&2; }
md5of(){ if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | awk '{print $1}'; else md5 -q "$1"; fi; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir)       OUTDIR="$2"; shift 2;;
    --release)      RELEASE="$2"; shift 2;;
    --build)        BUILD="$2"; shift 2;;
    --species)      SPECIES="$2"; shift 2;;
    --record-local) RECORD_LOCAL="$2"; shift 2;;
    -h|--help)      grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
    gtf|regulatory|rmsk|simplerepeat|segdups|dbsnp|genome|all) TARGETS+=("$1"); shift;;
    *) err "unknown argument: $1";;
  esac
done
[[ ${#TARGETS[@]} -eq 0 && -z "$RECORD_LOCAL" ]] && TARGETS=("gtf")
[[ " ${TARGETS[*]:-} " == *" all "* ]] && TARGETS=(genome gtf regulatory rmsk simplerepeat segdups dbsnp)

mkdir -p "$OUTDIR"
PROV="$OUTDIR/refs.provenance.json"
JSONL="$OUTDIR/refs.provenance.jsonl"
TS="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
ENSEMBL="https://ftp.ensembl.org/pub/release-${RELEASE}"
UCSC="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database"

# Seed the append-only JSONL from any pre-existing JSON manifest (one-time migration).
if [[ -f "$PROV" && ! -f "$JSONL" ]] && command -v python3 >/dev/null 2>&1; then
  python3 - "$PROV" "$JSONL" <<'PY'
import json,sys
try:
    d=json.load(open(sys.argv[1]))
    with open(sys.argv[2],"w") as o:
        for e in d.get("references",[]): o.write(json.dumps(e)+"\n")
except Exception: pass
PY
fi

record() {  # name url localpath
  local name="$1" url="$2" path="$3" md5 bytes
  md5="$(md5of "$path")"
  bytes="$(stat -f%z "$path" 2>/dev/null || stat -c%s "$path")"
  printf '{"name":"%s","build":"%s","release":"%s","source":"%s","path":"%s","bytes":%s,"md5":"%s","fetched":"%s"}\n' \
    "$name" "$BUILD" "$RELEASE" "$url" "$path" "$bytes" "$md5" "$TS" >> "$JSONL"
  log "✓ $name  md5=$md5  ($bytes bytes)"
}

fetch() {  # url localpath
  local url="$1" path="$2"
  [[ -s "$path" ]] && { log "exists, skip: $path"; return; }
  log "downloading: $url"
  curl -fsSL --retry 3 -o "$path" "$url" || err "download failed: $url"
}

if [[ -n "$RECORD_LOCAL" ]]; then
  [[ -s "$RECORD_LOCAL" ]] || err "no such file: $RECORD_LOCAL"
  log "checksumming local reference (whole genome may take a moment)…"
  record "genome_fasta_local" "file://$RECORD_LOCAL" "$RECORD_LOCAL"
fi

for t in "${TARGETS[@]:-}"; do case "$t" in
  genome)
    f="$OUTDIR/Homo_sapiens.${BUILD}.dna.primary_assembly.fa.gz"
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
  simplerepeat)
    f="$OUTDIR/simpleRepeat.${BUILD}.txt.gz"
    fetch "${UCSC}/simpleRepeat.txt.gz" "$f"
    record "simplerepeat" "${UCSC}/simpleRepeat.txt.gz" "$f" ;;
  segdups)
    f="$OUTDIR/genomicSuperDups.${BUILD}.txt.gz"
    fetch "${UCSC}/genomicSuperDups.txt.gz" "$f"
    record "segdups" "${UCSC}/genomicSuperDups.txt.gz" "$f" ;;
  dbsnp)
    f="$OUTDIR/${SPECIES}.dbsnp.${BUILD}.vcf.gz"
    fetch "${ENSEMBL}/variation/vcf/${SPECIES}/${SPECIES}.vcf.gz" "$f"
    record "dbsnp_vcf" "${ENSEMBL}/variation/vcf/${SPECIES}/" "$f" ;;
esac; done

# ── rebuild merged provenance JSON (dedup by name, latest wins) ─────────────
if command -v python3 >/dev/null 2>&1 && [[ -f "$JSONL" ]]; then
  python3 - "$JSONL" "$PROV" "$BUILD" "$RELEASE" "$TS" <<'PY'
import json,sys
jsonl,out,build,release,ts=sys.argv[1:6]
ents={}
for line in open(jsonl):
    line=line.strip()
    if not line: continue
    try: e=json.loads(line); ents[e["name"]]=e
    except Exception: pass
json.dump({"generated":ts,"build":build,"ensembl_release":release,
           "references":list(ents.values())}, open(out,"w"), indent=2)
PY
  log "provenance manifest updated → $PROV"
fi
log "DONE."
