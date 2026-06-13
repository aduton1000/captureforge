#!/usr/bin/env python3
"""
CaptureForge — test fixture : make_test_reference.py

Build a tiny self-contained mini-genome from the full reference for the `test`
profile: each test gene's locus (± flank) becomes its own contig `<chrom>_<gene>`
with LOCAL coordinates, and the GTF / rsID anchors are coordinate-shifted to
match. This lets the full 10-stage pipeline run end-to-end on a Mac in minutes
without subsetting whole (huge) chromosomes.

Usage:
  make_test_reference.py --genes test_genes.csv --gtf full.gtf.gz --fasta full.fa \
     [--rsid-bed full_anchors.bed] --flank 30000 --outdir test/data/ref
"""
from __future__ import annotations
import argparse, csv, gzip, os, re, subprocess, sys
from pathlib import Path

NAME_RE = re.compile(r'gene_name "([^"]+)"')


def read_genes(path):
    genes = []
    if path.endswith((".csv", ".tsv")):
        d = "\t" if path.endswith(".tsv") else ","
        for r in csv.DictReader(open(path), delimiter=d):
            genes.append(r["gene"].strip())
    else:
        for l in open(path):
            s = l.strip().split(",")[0].split("\t")[0]
            if s and not s.startswith("#") and s.lower() != "gene":
                genes.append(s)
    return genes


def gene_spans(gtf, wanted):
    spans = {}
    op = gzip.open if gtf.endswith(".gz") else open
    with op(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            m = NAME_RE.search(f[8])
            if m and m.group(1) in wanted:
                spans[m.group(1)] = (f[0], int(f[3]), int(f[4]))
    return spans


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--rsid-bed", default=None)
    ap.add_argument("--flank", type=int, default=30000)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    od = Path(args.outdir); od.mkdir(parents=True, exist_ok=True)
    genes = read_genes(args.genes)
    spans = gene_spans(args.gtf, set(genes))
    missing = [g for g in genes if g not in spans]
    if missing:
        sys.exit(f"ERROR: genes not found in GTF: {missing}")

    # per-gene region + new contig name + offset (region_start, 1-based)
    region = {}   # gene -> (chrom, rstart, rend, contig)
    for g in genes:
        ch, s, e = spans[g]
        rs = max(1, s - args.flank); re_ = e + args.flank
        region[g] = (ch, rs, re_, f"{ch}_{g}")

    # ── mini FASTA via samtools faidx ──
    fa = od / "test_genome.fa"
    with open(fa, "w") as out:
        for g in genes:
            ch, rs, re_, contig = region[g]
            seq = subprocess.run(["samtools", "faidx", args.fasta, f"{ch}:{rs}-{re_}"],
                                 capture_output=True, text=True, check=True).stdout
            body = "".join(seq.splitlines()[1:])
            out.write(f">{contig}\n")
            for i in range(0, len(body), 60):
                out.write(body[i:i+60] + "\n")
    subprocess.run(["samtools", "faidx", str(fa)], check=True)

    # ── mini GTF (shift coords; rename seqname) ──
    op = gzip.open if args.gtf.endswith(".gz") else open
    gtf_out = od / "test_genes.gtf"
    # map chrom -> list of (rstart, rend, contig) for membership test
    by_chrom = {}
    for g, (ch, rs, re_, contig) in region.items():
        by_chrom.setdefault(ch, []).append((rs, re_, contig, g))
    wanted = set(genes)
    with op(args.gtf, "rt") as fh, open(gtf_out, "w") as out:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            m = NAME_RE.search(f[8])
            if not m or m.group(1) not in wanted:
                continue
            g = m.group(1); ch, rs, re_, contig = region[g]
            s, e = int(f[3]), int(f[4])
            if f[0] != ch or e < rs or s > re_:
                continue
            f[0] = contig; f[3] = str(s - rs + 1); f[4] = str(e - rs + 1)
            out.write("\t".join(f) + "\n")

    # ── tiny regulatory GFF (one enhancer overlapping each gene's 5' region) ──
    # lets the `test` profile exercise the regulatory design class without the
    # full Ensembl Regulatory Build.
    reg = od / "test_regulatory.gff"
    with open(reg, "w") as out:
        for g in genes:
            ch, rs, re_, contig = region[g]
            gs_mini = spans[g][1] - rs + 1
            out.write(f"{contig}\tCaptureForgeTest\tenhancer\t{gs_mini+100}\t{gs_mini+700}\t"
                      f".\t+\t.\tID=enh_{g}\n")

    # ── shifted rsID anchors ──
    if args.rsid_bed and os.path.exists(args.rsid_bed):
        rb = od / "test_anchors.rsid.bed"
        with open(args.rsid_bed) as fh, open(rb, "w") as out:
            for line in fh:
                f = line.split()
                if len(f) < 4:
                    continue
                ch, pos0 = f[0], int(f[1])
                for rs, re_, contig, g in by_chrom.get(ch, []):
                    if rs <= pos0 + 1 <= re_:
                        out.write(f"{contig}\t{pos0-rs+1}\t{pos0-rs+2}\t{f[3]}\n")
                        break

    # ── BLAST db for specificity (no-alt mini-genome) ──
    subprocess.run(["makeblastdb", "-in", str(fa), "-dbtype", "nucl",
                    "-out", str(od / "test_noalt"), "-title", "test_noalt"],
                   check=True, capture_output=True)

    sys.stderr.write(f"[make_test_reference] {len(genes)} genes → {fa} "
                     f"({fa.stat().st_size//1024} KB), GTF + anchors + BLAST db in {od}\n")


if __name__ == "__main__":
    main()
