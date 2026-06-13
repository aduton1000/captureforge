#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════════
#  CaptureForge — bin/build_manual.sh
#  Render docs/manual/captureforge_manual.md → .pdf (xelatex) and .docx (pandoc).
#  Run from anywhere; builds inside docs/manual/ so figures/ paths resolve.
#  PDF engine: xelatex (preferred). Fallback: weasyprint (HTML→PDF) if present.
# ═══════════════════════════════════════════════════════════════════════════
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DIR="$ROOT/docs/manual"
MD="captureforge_manual.md"
cd "$DIR"

command -v pandoc >/dev/null || { echo "ERROR: pandoc not found" >&2; exit 1; }

COMMON=(--toc --number-sections --resource-path=.:figures)

echo "[build_manual] DOCX …"
pandoc "$MD" -o captureforge_manual.docx "${COMMON[@]}"

echo "[build_manual] PDF (xelatex) …"
if command -v xelatex >/dev/null 2>&1; then
  pandoc "$MD" -o captureforge_manual.pdf --pdf-engine=xelatex "${COMMON[@]}" \
    --include-in-header pdf-header.tex \
    -V geometry:margin=2.3cm -V colorlinks=true -V linkcolor=RoyalBlue -V toccolor=black \
    -V mainfont="Helvetica Neue" -V monofont="Menlo" \
    || pandoc "$MD" -o captureforge_manual.pdf --pdf-engine=xelatex "${COMMON[@]}" \
         --include-in-header pdf-header.tex \
         -V geometry:margin=2.3cm -V colorlinks=true   # retry with default font
elif command -v weasyprint >/dev/null 2>&1; then
  echo "[build_manual] xelatex absent — falling back to weasyprint (HTML→PDF)"
  pandoc "$MD" -o /tmp/_cfman.html --self-contained "${COMMON[@]}" --metadata title="CaptureForge Manual"
  weasyprint /tmp/_cfman.html captureforge_manual.pdf
else
  echo "ERROR: no PDF engine (install TeX Live/TinyTeX for xelatex, or weasyprint)" >&2; exit 1
fi

echo "[build_manual] done:"
ls -lh captureforge_manual.pdf captureforge_manual.docx | sed 's/^/  /'
