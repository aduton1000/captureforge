#!/usr/bin/env python3
"""
CaptureForge — Stage 7 : boost_dedup.py

1. Boosting: add extra COPIES over boost-flagged regions (promoter/regulatory,
   anchors, high-GC promoters). Copies increase a probe's molar representation in
   the pool; they do NOT add distinct oligos.
2. Dedup: collapse identical bait SEQUENCES across overlapping targets, summing
   copy numbers and recording all source regions.

Pool size is the number of UNIQUE sequences (an array-synth pool has one feature
per distinct oligo). The cost / pool-tier figure is therefore unique-based.

Tripwire: if unique baits exceed the smallest configured pool tier, the summary
sets `tripwire_exceeds_tier` and the script warns loudly — do NOT silently spill
into a larger tier.

Outputs: baits.dedup.tsv, baits.dedup.bed, baits.dedup.fasta, stage7_summary.json
"""
from __future__ import annotations
import argparse, json, re, sys
from collections import defaultdict
from pathlib import Path


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = re.sub(r"\([+-]\)$", "", line[1:].strip().split()[0])
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def gc_frac(s):
    return (s.count("G") + s.count("C")) / len(s) if s else 0.0


def base_copies(cls, gc, args):
    """Boost matrix: promoter/regulatory and anchor get +1; high-GC promoter +1 more."""
    n = 1
    if cls in ("promoter", "regulatory"):
        n += 1
        if gc >= args.high_gc:
            n += 1
    if cls == "anchor":
        n += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", required=True, help="retiled passing baits BED")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--high-gc", type=float, default=0.65)
    ap.add_argument("--per-oligo-price", type=float, default=0.00200)
    ap.add_argument("--pool-tiers", default="12000,92000,244000")
    ap.add_argument("--prefix", default="baits.dedup")
    ap.add_argument("--summary", default=None)
    args = ap.parse_args()

    seqs = read_fasta(args.fasta)
    tiers = [int(x) for x in args.pool_tiers.split(",")]

    # group by identical sequence (dedup); join bed<->fasta on stripped bait_id
    groups = {}            # seq -> dict(copies, sources[], classes set, genes set, rep_line)
    n_in = 0
    missing = 0
    for line in open(args.bed):
        if not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        chrom, s, e, name = f[0], f[1], f[2], f[3]
        strand = f[5] if len(f) > 5 else "+"
        parts = name.split("|")
        gene, cls = parts[0], (parts[1] if len(parts) > 1 else "coding")
        seq = seqs.get(name)
        if not seq:
            missing += 1
            continue
        n_in += 1
        gc = gc_frac(seq)
        copies = base_copies(cls, gc, args)
        g = groups.setdefault(seq, {"copies": 0, "sources": [], "classes": set(),
                                    "genes": set(), "rep": (chrom, s, e, strand), "gc": gc})
        g["copies"] += copies
        g["sources"].append(f"{chrom}:{s}-{e}:{name}")
        g["classes"].add(cls); g["genes"].add(gene)

    # write outputs
    tsv = open(f"{args.prefix}.tsv", "w")
    bed = open(f"{args.prefix}.bed", "w")
    fa = open(f"{args.prefix}.fasta", "w")
    tsv.write("bait_id\tsequence\tcopy_number\tn_sources\tgenes\tclasses\tgc\tsource_regions\n")
    copydist = defaultdict(int)
    i = 0
    for seq, g in groups.items():
        i += 1
        chrom, s, e, strand = g["rep"]
        genes = ",".join(sorted(g["genes"])); classes = ",".join(sorted(g["classes"]))
        bait_id = f"CF_{i:06d}|{sorted(g['genes'])[0]}|{sorted(g['classes'])[0]}"
        tsv.write(f"{bait_id}\t{seq}\t{g['copies']}\t{len(g['sources'])}\t{genes}\t{classes}\t"
                  f"{g['gc']:.3f}\t{';'.join(g['sources'])}\n")
        bed.write(f"{chrom}\t{s}\t{e}\t{bait_id}\t{g['copies']}\t{strand}\n")
        fa.write(f">{bait_id} copies={g['copies']} genes={genes} classes={classes}\n{seq}\n")
        copydist[g["copies"]] += 1
    tsv.close(); bed.close(); fa.close()

    unique = len(groups)
    total_oligo_molecules = sum(g["copies"] for g in groups.values())
    fit = next((t for t in tiers if unique <= t), None)
    cost = unique * args.per_oligo_price

    summary = {
        "n_baits_in": n_in,
        "n_unique": unique,
        "total_copies": total_oligo_molecules,
        "copy_number_distribution": {str(k): v for k, v in sorted(copydist.items())},
        "synthesis_pool_cost_usd": round(cost, 2),
        "per_oligo_price": args.per_oligo_price,
        "pool_tiers": tiers,
        "pool_tier_fit": fit,
        "tripwire_exceeds_tier": unique > tiers[0],
        "missing_sequence": missing,
    }
    if args.summary:
        Path(args.summary).write_text(json.dumps(summary, indent=2))
    sys.stderr.write(f"[boost_dedup] unique={unique} (in={n_in}) total_copies={total_oligo_molecules} "
                     f"tier_fit={fit} cost=${cost:.2f}\n")
    if summary["tripwire_exceeds_tier"]:
        sys.stderr.write(f"[boost_dedup] ⚠ TRIPWIRE: unique baits {unique} EXCEED smallest pool tier "
                         f"{tiers[0]} — do not spill silently; decide trim-vs-pay.\n")


if __name__ == "__main__":
    main()
