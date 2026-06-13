#!/usr/bin/env python3
"""
CaptureForge — Stage 6 : gap_retile.py

Find masked-target regions not covered by passing baits, then iteratively try to
close them: tile denser candidates over each gap, re-check specificity, and
re-check composition with RELAXED bounds (within documented limits). Iterate up
to --max-passes or until the coverage target is met. Residual gaps are reported
WITH a reason each (never silently dropped).

Reuses the verified Stage 4/5 tools as subprocesses (specificity_filter.py,
composition_filter.py) so the same logic governs re-tiled baits.

Outputs: retiled.pass.bed, retiled.pass.fasta, stage6_summary.json,
residual_gaps.bed (chrom start end gene|class reason).
"""
from __future__ import annotations
import argparse, json, math, os, shutil, subprocess, sys, tempfile
from collections import defaultdict
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))


def tool(name):
    """Resolve a sibling CLI: prefer co-located (local), else PATH (Nextflow stages bin/)."""
    p = os.path.join(HERE, name)
    return p if os.path.exists(p) else (shutil.which(name) or name)


def merge(ivs):
    ivs = sorted(ivs); out = []
    for s, e in ivs:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def load_class_bed(path):
    d = defaultdict(list)               # class -> list[(chrom,s,e,gene,name)]
    for l in open(path):
        if not l.strip():
            continue
        f = l.rstrip("\n").split("\t")
        name = f[3]; parts = name.split("|")
        d[parts[1] if len(parts) > 1 else "coding"].append(
            (f[0], int(f[1]), int(f[2]), parts[0], name))
    return d


def footprint(bed_paths):
    fp = defaultdict(list)
    for p in bed_paths:
        for l in open(p):
            if l.strip():
                f = l.split("\t"); fp[f[0]].append((int(f[1]), int(f[2])))
    return {c: merge(v) for c, v in fp.items()}


def gaps_for(intervals, fpc):
    """intervals on one chrom minus merged footprint fpc -> gap subintervals."""
    fpc = merge(fpc); out = []
    for s, e in intervals:
        cur = [(s, e)]
        for fs, fe in fpc:
            nxt = []
            for cs, ce in cur:
                if fe <= cs or fs >= ce:
                    nxt.append((cs, ce)); continue
                if fs > cs:
                    nxt.append((cs, fs))
                if fe < ce:
                    nxt.append((fe, ce))
            cur = nxt
        out += [(cs, ce) for cs, ce in cur if ce > cs]
    return out


def coverage(masked_by_class, fp):
    cov = {}
    allm = allc = 0
    for cls, items in masked_by_class.items():
        bychrom = defaultdict(list)
        for ch, s, e, *_ in items:
            bychrom[ch].append((s, e))
        m = c = 0
        for ch, ivs in bychrom.items():
            ivs = merge(ivs)
            m += sum(e - s for s, e in ivs)
            g = gaps_for(ivs, fp.get(ch, []))
            c += sum(e - s for s, e in ivs) - sum(e - s for s, e in g)
        cov[cls] = {"masked_bp": m, "covered_bp": c, "cov": round(c / m, 4) if m else 1.0}
        allm += m; allc += c
    cov["OVERALL"] = {"masked_bp": allm, "covered_bp": allc,
                      "cov": round(allc / allm, 4) if allm else 1.0}
    return cov


def tile_gap(s, e, L, step, chrom_len):
    if e - s <= L:
        mid = (s + e) // 2
        bs = max(0, min(mid - L // 2, chrom_len - L))
        return [(bs, bs + L)]
    starts, bs = [], s
    while bs + L < e:
        starts.append(bs); bs += step
    starts.append(e - L)
    seen, out = set(), []
    for bs in starts:
        bs = max(0, min(bs, chrom_len - L))
        if bs not in seen:
            seen.add(bs); out.append((bs, bs + L))
    return out


def run(cmd):
    p = subprocess.run(cmd, text=True, capture_output=True)
    if p.returncode != 0:
        sys.stderr.write(f"[gap_retile] subprocess failed: {' '.join(map(str, cmd))}\n{p.stderr[-1500:]}\n")
        raise subprocess.CalledProcessError(p.returncode, cmd)
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--masked", required=True)
    ap.add_argument("--passing-bed", required=True)
    ap.add_argument("--passing-fasta", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--db", required=True)
    ap.add_argument("--aligner", default="blast")
    ap.add_argument("--bait-len", type=int, default=120)
    ap.add_argument("--max-passes", type=int, default=3)
    ap.add_argument("--coverage-target", type=float, default=0.95)
    ap.add_argument("--threads", type=int, default=4)
    # specificity thresholds (passed through)
    ap.add_argument("--min-identity", type=float, default=0.90)
    ap.add_argument("--min-lenfrac", type=float, default=0.50)
    ap.add_argument("--max-offtarget", type=int, default=1)
    # composition (relaxed per pass within these documented bounds)
    ap.add_argument("--gc-min", type=float, default=0.25)
    ap.add_argument("--gc-max", type=float, default=0.70)
    ap.add_argument("--gc-relax", type=float, default=0.05)
    ap.add_argument("--tm-min", type=float, default=60.0)
    ap.add_argument("--tm-max", type=float, default=90.0)
    ap.add_argument("--tm-relax", type=float, default=5.0)
    ap.add_argument("--homopolymer-max", type=int, default=10)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--fasta-out", required=True)
    ap.add_argument("--residual", required=True)
    ap.add_argument("--summary", default=None)
    args = ap.parse_args()

    sizes = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(args.fai)}
    masked = load_class_bed(args.masked)

    # seed working set from the Stage-5 passing baits
    work = tempfile.mkdtemp(prefix="retile_")
    cur_bed = Path(work) / "pass.bed"; cur_fa = Path(work) / "pass.fasta"
    cur_bed.write_bytes(Path(args.passing_bed).read_bytes())
    cur_fa.write_bytes(Path(args.passing_fasta).read_bytes())

    passes = []
    added_total = 0
    last_gap_verdicts = {}

    for p in range(1, args.max_passes + 1):
        fp = footprint([str(cur_bed)])
        cov = coverage(masked, fp)
        passes.append({"pass": p - 1, "coverage": {k: v["cov"] for k, v in cov.items()}})
        sys.stderr.write(f"[gap_retile] before pass {p}: overall cov={cov['OVERALL']['cov']:.3f}\n")
        if cov["OVERALL"]["cov"] >= args.coverage_target:
            break

        # build gap candidates (denser step = L/3)
        step = max(1, round(args.bait_len / 3))
        gapcand = Path(work) / f"gapcand_p{p}.bed"
        n = defaultdict(int)
        with open(gapcand, "w") as o:
            for cls, items in masked.items():
                bychrom = defaultdict(list)
                meta = {}
                for ch, s, e, gene, name in items:
                    bychrom[ch].append((s, e)); meta[(ch, s, e)] = gene
                for ch, ivs in bychrom.items():
                    for gs, ge in gaps_for(merge(ivs), fp.get(ch, [])):
                        gene = next((g for (c2, s2, e2), g in meta.items()
                                     if c2 == ch and s2 <= gs and e2 >= ge), "NA")
                        for bs, be in tile_gap(gs, ge, args.bait_len, step, sizes.get(ch, ge + args.bait_len)):
                            n[(gene, cls)] += 1
                            o.write(f"{ch}\t{bs}\t{be}\t{gene}|{cls}|gapret|p{p}b{n[(gene,cls)]:04d}\t0\t+\n")
        if gapcand.stat().st_size == 0:
            break
        gapfa = Path(work) / f"gapcand_p{p}.fasta"
        run(["bedtools", "getfasta", "-fi", args.fasta, "-bed", str(gapcand), "-nameOnly", "-fo", str(gapfa)])

        # specificity
        spec_bed = Path(work) / f"gapspec_p{p}.bed"; spec_tsv = Path(work) / f"gapspec_p{p}.tsv"
        run([sys.executable, tool("specificity_filter.py"),
             "--candidates", str(gapfa), "--candidates-bed", str(gapcand),
             "--aligner", args.aligner, "--db", args.db,
             "--min-identity", str(args.min_identity), "--min-lenfrac", str(args.min_lenfrac),
             "--max-offtarget", str(args.max_offtarget), "--threads", str(args.threads),
             "--no-ref-check",   # main SPECIFICITY_FILTER already validated the reference
             "-o", str(spec_bed), "--verdicts", str(spec_tsv)])

        # composition with RELAXED bounds
        comp_pre = Path(work) / f"gapcomp_p{p}"
        run([sys.executable, tool("composition_filter.py"),
             "--pass-bed", str(spec_bed), "--fasta", str(gapfa),
             "--gc-min", str(args.gc_min - args.gc_relax), "--gc-max", str(args.gc_max + args.gc_relax),
             "--gc-max-promoter", str(args.gc_max + args.gc_relax + 0.10),
             "--homopolymer-max", str(args.homopolymer_max + 1),
             "--tm-min", str(args.tm_min - args.tm_relax), "--tm-max", str(args.tm_max + args.tm_relax),
             "--prefix", str(comp_pre), "--verdicts", str(comp_pre) + ".tsv"])

        # record verdicts for residual-gap reasoning
        last_gap_verdicts = {"spec": str(spec_tsv), "comp": str(comp_pre) + ".tsv"}

        add_bed = Path(str(comp_pre) + ".pass.bed"); add_fa = Path(str(comp_pre) + ".pass.fasta")
        n_add = sum(1 for _ in open(add_bed)) if add_bed.exists() else 0
        if n_add == 0:
            sys.stderr.write(f"[gap_retile] pass {p}: no new baits recovered → stop\n")
            break
        with open(cur_bed, "a") as o:
            o.write(add_bed.read_text())
        with open(cur_fa, "a") as o:
            o.write(add_fa.read_text())
        added_total += n_add
        sys.stderr.write(f"[gap_retile] pass {p}: +{n_add} baits recovered\n")

    # finalize
    fp = footprint([str(cur_bed)])
    final_cov = coverage(masked, fp)
    Path(args.output).write_bytes(cur_bed.read_bytes())
    Path(args.fasta_out).write_bytes(cur_fa.read_bytes())

    # residual gaps + reason
    reasons = defaultdict(int)
    with open(args.residual, "w") as o:
        for cls, items in masked.items():
            bychrom = defaultdict(list); meta = {}
            for ch, s, e, gene, name in items:
                bychrom[ch].append((s, e)); meta.setdefault(ch, []).append((s, e, gene))
            for ch, ivs in bychrom.items():
                for gs, ge in gaps_for(merge(ivs), fp.get(ch, [])):
                    gene = next((g for (s2, e2, g) in meta[ch] if s2 <= gs and e2 >= ge), "NA")
                    w = ge - gs
                    reason = ("edge_short" if w < args.bait_len * 0.3
                              else "non_specific_or_composition")
                    reasons[reason] += 1
                    o.write(f"{ch}\t{gs}\t{ge}\t{gene}|{cls}\t{reason}\n")

    converged = final_cov["OVERALL"]["cov"] >= args.coverage_target
    summary = {
        "passes": passes,
        "final_coverage": {k: v["cov"] for k, v in final_cov.items()},
        "final_coverage_detail": final_cov,
        "baits_recovered": added_total,
        "converged": converged,
        "coverage_target": args.coverage_target,
        "residual_gap_reasons": dict(reasons),
        "n_baits_final": sum(1 for _ in open(args.output)),
    }
    if args.summary:
        Path(args.summary).write_text(json.dumps(summary, indent=2))
    sys.stderr.write(f"[gap_retile] final overall cov={final_cov['OVERALL']['cov']:.3f} "
                     f"converged={converged} recovered={added_total}\n")


if __name__ == "__main__":
    main()
