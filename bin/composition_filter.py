#!/usr/bin/env python3
"""
CaptureForge — Stage 5 : composition_filter.py

Filter baits on sequence composition: GC%, homopolymer run, low-complexity
(Shannon entropy), nearest-neighbour Tm (primer3), and optional ViennaRNA
self-structure ΔG. Promoter-class baits use a relaxed upper GC bound.

Outputs: composition.tsv (per-bait metrics + verdict), <prefix>.pass.bed,
<prefix>.pass.fasta, summary JSON (histograms + rejection-by-reason). Also
reports an early unique-bait count (identical sequences collapsed).
"""
from __future__ import annotations
import argparse, json, math, re, sys
from collections import defaultdict, Counter
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
    n = len(s)
    return (s.count("G") + s.count("C")) / n if n else 0.0


def max_homopolymer(s):
    best = cur = 1
    for i in range(1, len(s)):
        cur = cur + 1 if s[i] == s[i - 1] else 1
        best = max(best, cur)
    return best if s else 0


def shannon_entropy(s):
    if not s:
        return 0.0
    c = Counter(s)
    n = len(s)
    return -sum((v / n) * math.log2(v / n) for v in c.values())


def tm_calc(s):
    try:
        import primer3
        return float(primer3.calc_tm(s))
    except Exception:
        try:
            from Bio.SeqUtils import MeltingTemp as mt
            return float(mt.Tm_NN(s))
        except Exception:
            return float("nan")


def dg_calc(s, enabled):
    if not enabled:
        return None
    try:
        import RNA
        _, mfe = RNA.fold(s.replace("T", "U"))
        return float(mfe)
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pass-bed", required=True, help="specificity-passing candidates.bed")
    ap.add_argument("--fasta", required=True, help="candidates.fasta (all baits)")
    ap.add_argument("--gc-min", type=float, default=0.25)
    ap.add_argument("--gc-max", type=float, default=0.70)
    ap.add_argument("--gc-max-promoter", type=float, default=0.80)
    ap.add_argument("--homopolymer-max", type=int, default=10)
    ap.add_argument("--tm-min", type=float, default=60.0)
    ap.add_argument("--tm-max", type=float, default=90.0)
    ap.add_argument("--min-entropy", type=float, default=1.2)
    ap.add_argument("--use-vienna", action="store_true")
    ap.add_argument("--deltaG-min", type=float, default=-9.0)
    ap.add_argument("--prefix", default="composition")
    ap.add_argument("--verdicts", default="composition.tsv")
    ap.add_argument("--summary", default=None)
    args = ap.parse_args()

    seqs = read_fasta(args.fasta)
    bed_lines = [l for l in open(args.pass_bed) if l.strip()]

    vout = open(args.verdicts, "w")
    vout.write("bait_id\tclass\tgc\ttm\thomopolymer\tentropy\tdeltaG\tverdict\treason\n")
    pass_bed = open(f"{args.prefix}.pass.bed", "w")
    pass_fa = open(f"{args.prefix}.pass.fasta", "w")

    by_class_total = defaultdict(int)
    by_class_reject = defaultdict(int)
    reasons = defaultdict(int)
    gc_pre, gc_post, tm_post = [], [], []
    uniq = {}
    n_pass = 0

    for line in bed_lines:
        f = line.rstrip("\n").split("\t")
        bait_id = f[3]
        cls = bait_id.split("|")[1] if "|" in bait_id else "coding"
        seq = seqs.get(bait_id, "")
        by_class_total[cls] += 1
        if not seq:
            continue
        gc = gc_frac(seq); hp = max_homopolymer(seq); ent = shannon_entropy(seq)
        tm = tm_calc(seq); dg = dg_calc(seq, args.use_vienna)
        gc_pre.append(gc)

        gc_hi = args.gc_max_promoter if cls in ("promoter", "regulatory") else args.gc_max
        why = []
        if not (args.gc_min <= gc <= gc_hi):
            why.append("GC")
        if hp > args.homopolymer_max:
            why.append("homopolymer")
        if ent < args.min_entropy:
            why.append("low_complexity")
        if math.isnan(tm) or not (args.tm_min <= tm <= args.tm_max):
            why.append("Tm")            # NaN Tm (e.g. N-containing bait) can't be validated → reject
        if dg is not None and dg < args.deltaG_min:
            why.append("structure")

        verdict = "PASS" if not why else "REJECT"
        vout.write(f"{bait_id}\t{cls}\t{gc:.3f}\t{tm:.1f}\t{hp}\t{ent:.2f}\t"
                   f"{'' if dg is None else round(dg,1)}\t{verdict}\t{','.join(why)}\n")
        for r in why:
            reasons[r] += 1
        if verdict == "PASS":
            pass_bed.write(line)
            pass_fa.write(f">{bait_id}\n{seq}\n")
            gc_post.append(gc); tm_post.append(tm)
            uniq.setdefault(seq, bait_id)
            n_pass += 1
        else:
            by_class_reject[cls] += 1

    vout.close(); pass_bed.close(); pass_fa.close()

    def hist(vals, lo, hi, n=20):
        vals = [v for v in vals if v == v]      # drop NaN defensively
        if not vals:
            return {}
        w = (hi - lo) / n
        h = defaultdict(int)
        for v in vals:
            b = min(n - 1, max(0, int((v - lo) / w))) if w else 0
            h[round(lo + b * w, 3)] += 1
        return dict(sorted(h.items()))

    summary = {
        "n_in": len(bed_lines), "n_pass": n_pass,
        "n_reject": sum(by_class_reject.values()),
        "reject_by_class": dict(by_class_reject),
        "total_by_class": dict(by_class_total),
        "reject_by_reason": dict(reasons),
        "unique_baits_after_dedup": len(uniq),
        "gc_hist_pre": hist(gc_pre, 0, 1), "gc_hist_post": hist(gc_post, 0, 1),
        "tm_hist_post": hist(tm_post, 40, 100),
        "params": {"gc": [args.gc_min, args.gc_max], "gc_promoter_hi": args.gc_max_promoter,
                   "homopolymer_max": args.homopolymer_max, "tm": [args.tm_min, args.tm_max],
                   "min_entropy": args.min_entropy, "use_vienna": args.use_vienna},
    }
    if args.summary:
        Path(args.summary).write_text(json.dumps(summary, indent=2))
    sys.stderr.write(f"[composition] pass={n_pass}/{len(bed_lines)} | unique(after dedup)={len(uniq)} | "
                     f"reject_by_reason={dict(reasons)}\n")


if __name__ == "__main__":
    main()
