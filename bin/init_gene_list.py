#!/usr/bin/env python3
"""
CaptureForge — bin/init_gene_list.py

Scaffold an EDITABLE design-spec CSV from a bare list of gene symbols. This is
the recommended starting point for a new panel: run it, then open the CSV and
add promoter / regulatory / cnv / str / anchor intent per gene before running
the pipeline with `--input <that csv>`.

Usage:
    init_gene_list.py genes.txt -o my_panel.spec.csv \
        [--default-classes coding,promoter]

Optional (suggestions the user CONFIRMS, never silent):
    --gtf <gtf> --segdups <bed> --rmsk <bed>   pre-flag cnv/str candidates
    (requires bedtools on PATH; off by default)
"""
from __future__ import annotations
import argparse, csv, sys
from pathlib import Path

COLUMNS = ["gene", "promoter", "regulatory", "cnv", "str", "extra_regions", "anchors"]

def read_symbols(path: str) -> list[str]:
    out, seen = [], set()
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            # tolerate CSV/TSV first column too
            s = s.split(",")[0].split("\t")[0].strip()
            if s and not s.startswith("#") and s.lower() != "gene" and s not in seen:
                seen.add(s); out.append(s)
    return out

def main():
    ap = argparse.ArgumentParser(description="Scaffold an editable CaptureForge spec CSV from bare symbols.")
    ap.add_argument("symbols", help="plain-text file, one gene symbol per line")
    ap.add_argument("-o", "--output", default="panel.spec.csv")
    ap.add_argument("--default-classes", default="coding,promoter",
                    help="classes applied to every gene by default (coding always implied)")
    ap.add_argument("--gtf", help="(optional) GTF for track-based cnv/str suggestions")
    ap.add_argument("--segdups", help="(optional) segmental-dup BED → suggest cnv")
    ap.add_argument("--rmsk", help="(optional) tandem-repeat BED → suggest str")
    args = ap.parse_args()

    genes = read_symbols(args.symbols)
    if not genes:
        sys.exit(f"ERROR: no gene symbols found in {args.symbols}")

    defaults = {c.strip() for c in args.default_classes.split(",")}
    promoter_default = 1 if "promoter" in defaults else 0

    suggestions = {}
    if args.segdups or args.rmsk:
        suggestions = suggest_from_tracks(genes, args.gtf, args.segdups, args.rmsk)

    with open(args.output, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(COLUMNS)
        for g in genes:
            sug = suggestions.get(g, {})
            w.writerow([g, promoter_default, 0,
                        sug.get("cnv", 0), sug.get("str", 0), "", ""])

    sys.stderr.write(
        f"[init_gene_list] wrote {len(genes)} genes → {args.output}\n"
        f"  Edit it to add intent: promoter/regulatory/cnv/str (1/0), "
        f"extra_regions (chr:start-end;…), anchors (rsIDs).\n"
        f"  coding is always implied. Then run:\n"
        f"    nextflow run main.nf -profile mac_local,conda --input {args.output} \\\n"
        f"        --genome_fasta <ref.fa> --gtf <ann.gtf>\n")
    if suggestions:
        flagged = {g: s for g, s in suggestions.items() if s}
        sys.stderr.write(f"  NOTE: pre-flagged {len(flagged)} genes from tracks "
                         f"(SUGGESTIONS — review before running): {sorted(flagged)}\n")


def suggest_from_tracks(genes, gtf, segdups, rmsk):
    """Best-effort cnv/str pre-flags from segdup/repeat overlap. Suggestions only."""
    import shutil, subprocess, tempfile, re
    if not gtf:
        sys.stderr.write("  (track suggestions need --gtf for gene coords; skipping)\n")
        return {}
    if not shutil.which("bedtools"):
        sys.stderr.write("  (bedtools not found; skipping track suggestions)\n")
        return {}
    # build a gene BED from the GTF
    import gzip
    op = gzip.open if gtf.endswith(".gz") else open
    want = set(genes)
    gene_bed = []
    name_re = re.compile(r'gene_name "([^"]+)"')
    with op(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            m = name_re.search(f[8])
            if m and m.group(1) in want:
                gene_bed.append((f[0], int(f[3]) - 1, int(f[4]), m.group(1)))
    out = {}
    with tempfile.TemporaryDirectory() as td:
        gb = Path(td) / "genes.bed"
        gb.write_text("".join(f"{c}\t{s}\t{e}\t{n}\n" for c, s, e, n in sorted(gene_bed)))
        for track, flag in ((segdups, "cnv"), (rmsk, "str")):
            if not track:
                continue
            try:
                r = subprocess.run(["bedtools", "intersect", "-u", "-a", str(gb), "-b", track],
                                   capture_output=True, text=True, check=True)
                for ln in r.stdout.splitlines():
                    g = ln.split("\t")[3]
                    out.setdefault(g, {})[flag] = 1
            except subprocess.CalledProcessError as e:
                sys.stderr.write(f"  (bedtools intersect failed for {flag}: {e.stderr[:120]}; skipping)\n")
    return out


if __name__ == "__main__":
    main()
