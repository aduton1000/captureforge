#!/usr/bin/env python3
"""
CaptureForge — Stage 3 : tile_candidates.py

Tile candidate baits of fixed length across each masked target, applying the
per-design-class step / spacing policy, then extract their sequences.

Per class (from params.tiling; step = round(bait_len * step_frac)):
  coding      step bait_len/2   sliding, end-anchored
  promoter    step bait_len/3   denser
  regulatory  step bait_len/3
  cnv         step bait_len/2   EVEN spacing across the locus (uniform gaps)
  str         step bait_len/3   dense over the (post-mask) unique flanks
  anchor      step bait_len/2   guarantees >=2 baits across the rsID window

Sub-bait-length targets get one bait of bait_len centered on the target.

Outputs: candidates.bed, candidates.fasta, summary JSON.
"""
from __future__ import annotations
import argparse, json, math, subprocess, sys
from collections import defaultdict
from pathlib import Path

DEFAULT_TILING = {
    "coding":     {"step_frac": 0.50, "even": False},
    "promoter":   {"step_frac": 0.33, "even": False},
    "regulatory": {"step_frac": 0.33, "even": False},
    "cnv":        {"step_frac": 0.50, "even": True},
    "str":        {"step_frac": 0.33, "even": False},
    "anchor":     {"step_frac": 0.50, "even": False},
    "extra":      {"step_frac": 0.50, "even": False},
}


def load_fai(fai):
    return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai)}


def tile_interval(s, e, L, step, even, chrom_len):
    """Return list of (start,end) bait windows of length L covering [s,e]."""
    width = e - s
    if width <= L:                                   # center a single bait
        mid = (s + e) // 2
        bs = max(0, min(mid - L // 2, chrom_len - L))
        return [(bs, bs + L)]
    starts = []
    if even:
        n = max(2, math.ceil(width / step))          # n baits, uniform spacing
        last = e - L
        for i in range(n):
            bs = round(s + (last - s) * i / (n - 1))
            starts.append(bs)
    else:
        bs = s
        while bs + L < e:
            starts.append(bs)
            bs += step
        starts.append(e - L)                         # end-anchored final bait
    # clip + dedupe preserving order
    seen, out = set(), []
    for bs in starts:
        bs = max(0, min(bs, chrom_len - L))
        if bs not in seen:
            seen.add(bs); out.append((bs, bs + L))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--masked", required=True, help="targets.masked.bed")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--bait-len", type=int, default=120)
    ap.add_argument("--tiling-json", default=None, help="JSON of per-class {step_frac,even}")
    ap.add_argument("-o", "--output", required=True, help="candidates.bed")
    ap.add_argument("--fasta-out", required=True, help="candidates.fasta")
    ap.add_argument("--summary", default=None)
    args = ap.parse_args()

    L = args.bait_len
    tiling = dict(DEFAULT_TILING)
    if args.tiling_json:
        tiling.update(json.loads(args.tiling_json))
    sizes = load_fai(args.fai)

    counts_gene = defaultdict(int)
    counts_class = defaultdict(int)
    n_idx = defaultdict(int)
    lengths = set()

    with open(args.output, "w") as out:
        for line in open(args.masked):
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            chrom, s, e, name = f[0], int(f[1]), int(f[2]), f[3]
            strand = f[5] if len(f) > 5 else "+"
            parts = name.split("|")
            gene, cls = parts[0], (parts[1] if len(parts) > 1 else "coding")
            pol = tiling.get(cls, DEFAULT_TILING["coding"])
            step = max(1, round(L * pol["step_frac"]))
            for bs, be in tile_interval(s, e, L, step, pol["even"], sizes.get(chrom, e + L)):
                n_idx[(gene, cls)] += 1
                bait_id = f"{gene}|{cls}|{parts[2] if len(parts) > 2 else '.'}|b{n_idx[(gene,cls)]:04d}"
                out.write(f"{chrom}\t{bs}\t{be}\t{bait_id}\t0\t{strand}\n")
                counts_gene[gene] += 1
                counts_class[cls] += 1
                lengths.add(be - bs)

    # extract sequences (name from col4). Plus-strand (no -s) so FASTA headers
    # exactly match the BED names — downstream stages key on bait_id, and the
    # strand suffix bedtools adds with -s would silently break those lookups.
    # bait_strand handling / reverse-complement happens at Stage 8 per config.
    subprocess.run(["bedtools", "getfasta", "-fi", args.fasta, "-bed", args.output,
                    "-nameOnly", "-fo", args.fasta_out], check=True)

    total = sum(counts_class.values())
    summary = {
        "bait_len": L,
        "n_candidates": total,
        "by_class": dict(counts_class),
        "by_gene": dict(counts_gene),
        "length_uniformity": {"distinct_lengths": sorted(lengths),
                              "all_full_length": lengths == {L} if lengths else True},
    }
    if args.summary:
        Path(args.summary).write_text(json.dumps(summary, indent=2))
    sys.stderr.write(f"[tile_candidates] {total} candidate baits  by_class={dict(counts_class)}\n")


if __name__ == "__main__":
    main()
