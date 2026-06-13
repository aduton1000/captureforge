#!/usr/bin/env python3
"""
CaptureForge — Stage 9 : aggregate_qc.py

Aggregate all stage metrics, compute honest coverage, emit the analysis-pipeline
handoff, and run QC cross-checks.

KEY: final_covered_targets.bed is the AS-BUILT coverage — the intersection of the
raw targets with the footprint the final baits ACTUALLY cover. It is NOT the raw
target set. The downstream variant-analysis pipeline measures against this.

Coverage is reported two ways:
  mappable_cov = covered / masked-target  (design efficiency over tileable sequence)
  raw_cov      = covered / raw-target     (honest end-to-end; the deficit vs mappable
                 is sequence intentionally removed by masking — repeats/segdups/STR
                 cores/Alu — which is genomically uncapturable, by design)

For CNV genes, reports whether the mappable bins are evenly spaced enough to
support read-depth CNV calling (the decision-relevant metric, not raw_cov).

Outputs: final_covered_targets.bed, metrics.json, coverage_by_gene.tsv,
coverage_by_class.tsv, stage9_summary.json
"""
from __future__ import annotations
import argparse, csv, json, statistics, sys
from collections import defaultdict
from pathlib import Path


def merge(ivs):
    ivs = sorted(ivs); out = []
    for s, e in ivs:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def load_named(path):
    """class -> gene -> chrom -> list[(s,e)] ; plus chrom map for footprint."""
    d = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for l in open(path):
        if not l.strip():
            continue
        f = l.rstrip("\n").split("\t")
        p = f[3].split("|")
        cls = p[1] if len(p) > 1 else "coding"
        d[cls][p[0]][f[0]].append((int(f[1]), int(f[2])))
    return d


def footprint(path):
    fp = defaultdict(list)
    for l in open(path):
        if l.strip():
            f = l.split("\t"); fp[f[0]].append((int(f[1]), int(f[2])))
    return {c: merge(v) for c, v in fp.items()}


def intersect(ivs, fpc):
    fpc = merge(fpc); out = []
    for s, e in ivs:
        for fs, fe in fpc:
            lo = max(s, fs); hi = min(e, fe)
            if hi > lo:
                out.append((lo, hi))
    return merge(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--masked", required=True)
    ap.add_argument("--final-bed", required=True, help="baits.dedup.bed (final unique baits)")
    ap.add_argument("--summaries-dir", default=None, help="dir with stageN_summary.json")
    ap.add_argument("--per-oligo-price", type=float, default=0.00200)
    ap.add_argument("--pool-tiers", default="12000,92000,244000")
    ap.add_argument("--covered-vs-raw-min", type=float, default=0.90)
    ap.add_argument("--bait-len", type=int, default=120)
    ap.add_argument("--cnv-maxgap-bp", type=int, default=2000, help="max tolerated inter-bin gap for depth CNV")
    ap.add_argument("-o", "--covered-bed", required=True)
    ap.add_argument("--metrics", required=True)
    ap.add_argument("--by-gene", default="coverage_by_gene.tsv")
    ap.add_argument("--by-class", default="coverage_by_class.tsv")
    ap.add_argument("--summary", default=None)
    args = ap.parse_args()

    raw = load_named(args.raw)
    masked = load_named(args.masked)
    fp = footprint(args.final_bed)

    # ── AS-BUILT covered targets (intersection of raw targets with bait footprint) ──
    covered_out = open(args.covered_bed, "w")
    cov_by_cg = defaultdict(int)        # (cls,gene)->covered bp
    raw_by_cg = defaultdict(int)
    msk_by_cg = defaultdict(int)
    for cls, genes in raw.items():
        for gene, chroms in genes.items():
            for ch, ivs in chroms.items():
                ivs = merge(ivs)
                raw_by_cg[(cls, gene)] += sum(e - s for s, e in ivs)
                cov = intersect(ivs, fp.get(ch, []))
                cov_by_cg[(cls, gene)] += sum(e - s for s, e in cov)
                for s, e in cov:
                    covered_out.write(f"{ch}\t{s}\t{e}\t{gene}|{cls}\t0\t+\n")
    covered_out.close()
    covmsk_by_cg = defaultdict(int)     # (cls,gene)-> masked∩footprint bp
    for cls, genes in masked.items():
        for gene, chroms in genes.items():
            for ch, ivs in chroms.items():
                ivs = merge(ivs)
                msk_by_cg[(cls, gene)] += sum(e - s for s, e in ivs)
                cm = intersect(ivs, fp.get(ch, []))
                covmsk_by_cg[(cls, gene)] += sum(e - s for s, e in cm)

    # ── coverage tables ──
    classes = sorted({c for c, _ in raw_by_cg})
    by_class = {}
    for cls in classes:
        rb = sum(v for (c, g), v in raw_by_cg.items() if c == cls)
        cb = sum(v for (c, g), v in cov_by_cg.items() if c == cls)     # raw∩footprint
        mb = sum(v for (c, g), v in msk_by_cg.items() if c == cls)
        cmb = sum(v for (c, g), v in covmsk_by_cg.items() if c == cls)  # masked∩footprint
        by_class[cls] = {"raw_bp": rb, "masked_bp": mb, "covered_bp": cb, "covered_masked_bp": cmb,
                         "raw_cov": round(cb / rb, 4) if rb else 1.0,
                         "mappable_cov": round(cmb / mb, 4) if mb else 1.0}
    tot_raw = sum(v["raw_bp"] for v in by_class.values())
    tot_msk = sum(v["masked_bp"] for v in by_class.values())
    tot_cov = sum(v["covered_bp"] for v in by_class.values())
    tot_covm = sum(v["covered_masked_bp"] for v in by_class.values())
    overall = {"raw_bp": tot_raw, "masked_bp": tot_msk, "covered_bp": tot_cov,
               "covered_masked_bp": tot_covm,
               "raw_cov": round(tot_cov / tot_raw, 4) if tot_raw else 1.0,
               "mappable_cov": round(tot_covm / tot_msk, 4) if tot_msk else 1.0}

    with open(args.by_class, "w") as o:
        o.write("class\traw_bp\tmasked_bp\tcovered_bp\traw_cov\tmappable_cov\n")
        for cls in classes:
            v = by_class[cls]
            o.write(f"{cls}\t{v['raw_bp']}\t{v['masked_bp']}\t{v['covered_bp']}\t{v['raw_cov']}\t{v['mappable_cov']}\n")
    genes = sorted({g for _, g in raw_by_cg})
    gene_cov = {}
    with open(args.by_gene, "w") as o:
        o.write("gene\traw_bp\tcovered_bp\traw_cov\n")
        for g in genes:
            rb = sum(v for (c, gg), v in raw_by_cg.items() if gg == g)
            cb = sum(v for (c, gg), v in cov_by_cg.items() if gg == g)
            gene_cov[g] = round(cb / rb, 4) if rb else 1.0
            o.write(f"{g}\t{rb}\t{cb}\t{gene_cov[g]}\n")
    # fail-loud per-gene flags (never silently drop a gene)
    fully_uncovered = sorted(g for g, c in gene_cov.items() if c == 0.0)
    low_coverage = sorted(((g, c) for g, c in gene_cov.items() if 0.0 < c < 0.50), key=lambda x: x[1])
    if fully_uncovered:
        sys.stderr.write(f"[aggregate_qc] ⚠ FULLY-UNCOVERED GENES (0% — need decision): {fully_uncovered}\n")

    # ── CNV bin-spacing (read-depth callability) ──
    cnv_bins = defaultdict(list)
    for l in open(args.final_bed):
        f = l.split("\t"); p = f[3].split("|")
        if len(p) > 1 and p[1] == "cnv":
            cnv_bins[p[1] + "::" + l.split("\t")[3].split("|")[0]].append(int(f[1]))
    # regroup by gene from bait_id genes field (CF ids carry gene in field2)
    cnv_by_gene = defaultdict(list)
    for l in open(args.final_bed):
        f = l.rstrip("\n").split("\t"); p = f[3].split("|")
        # final bait_id = CF_xxxxxx|gene|class
        if len(p) >= 3 and p[2] == "cnv":
            cnv_by_gene[p[1]].append(int(f[1]))
    cnv_spacing = {}
    for gene, starts in cnv_by_gene.items():
        starts = sorted(starts)
        gaps = [b - a for a, b in zip(starts, starts[1:])]
        if not gaps:
            continue
        cnv_spacing[gene] = {
            "n_bins": len(starts),
            "median_gap": int(statistics.median(gaps)),
            "max_gap": max(gaps),
            "cv_gap": round(statistics.pstdev(gaps) / statistics.mean(gaps), 3) if statistics.mean(gaps) else 0,
            "depth_callable": max(gaps) <= args.cnv_maxgap_bp and len(starts) >= 10,
        }

    tiers = [int(x) for x in args.pool_tiers.split(",")]
    # pull unique count & cost from Stage 7 summary if present
    n_unique = None
    sdir = Path(args.summaries_dir) if args.summaries_dir else None
    stage_summaries = {}
    if sdir and sdir.exists():
        for jf in sorted(sdir.glob("stage*_summary.json")):
            try:
                stage_summaries[jf.stem] = json.loads(jf.read_text())
            except Exception:
                pass
        s7 = stage_summaries.get("stage7_summary", {})
        n_unique = s7.get("n_unique")
    cost = (n_unique or 0) * args.per_oligo_price
    fit = next((t for t in tiers if n_unique and n_unique <= t), None)

    raw_gate_pass = overall["raw_cov"] >= args.covered_vs_raw_min
    metrics = {
        "coverage_overall": overall,
        "coverage_by_class": by_class,
        "n_unique_baits": n_unique,
        "synthesis_pool_cost_usd": round(cost, 2),
        "pool_tier_fit": fit, "pool_tiers": tiers,
        "tripwire_exceeds_tier": bool(n_unique and n_unique > tiers[0]),
        "cnv_bin_spacing": cnv_spacing,
        "qc_gates": {
            "raw_cov_vs_min": {"value": overall["raw_cov"], "min": args.covered_vs_raw_min,
                               "pass": raw_gate_pass,
                               "note": ("raw_cov below threshold is EXPECTED — the deficit vs mappable_cov "
                                        f"({overall['mappable_cov']}) is masked/uncapturable sequence")},
            "mappable_cov": overall["mappable_cov"],
            "fully_uncovered_genes": fully_uncovered,
            "low_coverage_genes": [{"gene": g, "raw_cov": c} for g, c in low_coverage],
        },
        "stage_summaries": stage_summaries,
    }
    Path(args.metrics).write_text(json.dumps(metrics, indent=2))
    if args.summary:
        Path(args.summary).write_text(json.dumps({k: metrics[k] for k in
            ("coverage_overall", "coverage_by_class", "n_unique_baits", "synthesis_pool_cost_usd",
             "pool_tier_fit", "tripwire_exceeds_tier", "cnv_bin_spacing", "qc_gates")}, indent=2))
    sys.stderr.write(f"[aggregate_qc] mappable_cov={overall['mappable_cov']} raw_cov={overall['raw_cov']} "
                     f"unique={n_unique} cost=${cost:.2f} tier={fit}\n")
    sys.stderr.write(f"[aggregate_qc] final_covered_targets.bed = AS-BUILT covered target (handoff)\n")


if __name__ == "__main__":
    main()
