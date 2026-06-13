#!/usr/bin/env python3
"""
CaptureForge — visualization : make_figures.py

Generate all QC/design figures (PNG + SVG) into results/figures/. Reads the
stage BEDs and JSON summaries from a working directory. Figures:
  coverage_by_class, coverage_by_gene, gc_hist, tm_hist, offtarget_hist,
  copynumber_hist, convergence, per_class_summary, cost_pool_tier,
  coverage_tracks_page{N} (genome-browser-style per-gene panels, paginated).
"""
from __future__ import annotations
import argparse, csv, json, math, os
from collections import defaultdict, Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

CLASS_COLORS = {"coding": "#1f77b4", "promoter": "#ff7f0e", "regulatory": "#2ca02c",
                "cnv": "#d62728", "str": "#9467bd", "anchor": "#8c564b", "extra": "#7f7f7f"}
plt.rcParams.update({"figure.dpi": 110, "font.size": 9, "axes.grid": True,
                     "grid.alpha": 0.3, "axes.axisbelow": True})


def save(fig, figdir, name):
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(figdir, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)


def load_json(p):
    return json.loads(Path(p).read_text()) if os.path.exists(p) else {}


def load_stage(W, n):
    """Find stage-N summary across naming conventions (Nextflow vs dev)."""
    for cand in (f"{W}/summaries/stage{n}_summary.json", f"{W}/stage{n}_summary.json", f"{W}/stage{n}.json"):
        if os.path.exists(cand):
            return load_json(cand)
    return {}


def load_bed_by_gene(path, gene_field):
    d = defaultdict(list)
    if not os.path.exists(path):
        return d
    for l in open(path):
        if not l.strip():
            continue
        f = l.rstrip("\n").split("\t")
        parts = f[3].split("|")
        gene = parts[gene_field] if len(parts) > gene_field else parts[0]
        cls = parts[1] if len(parts) > 1 else "coding"
        d[gene].append((f[0], int(f[1]), int(f[2]), cls))
    return d


def fig_coverage_by_class(W, figdir):
    rows = list(csv.DictReader(open(os.path.join(W, "coverage_by_class.tsv")), delimiter="\t"))
    cls = [r["class"] for r in rows]
    raw = [float(r["raw_cov"]) * 100 for r in rows]
    mp = [float(r["mappable_cov"]) * 100 for r in rows]
    x = range(len(cls)); w = 0.38
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar([i - w/2 for i in x], mp, w, label="mappable cov", color="#2ca02c")
    ax.bar([i + w/2 for i in x], raw, w, label="raw cov", color="#1f77b4")
    ax.set_xticks(list(x)); ax.set_xticklabels(cls)
    ax.set_ylabel("% covered"); ax.set_ylim(0, 100)
    ax.axhline(95, ls="--", c="grey", lw=1)
    ax.set_title("Coverage by design class (raw vs mappable)"); ax.legend()
    save(fig, figdir, "coverage_by_class")


def fig_coverage_by_gene(W, figdir):
    rows = list(csv.DictReader(open(os.path.join(W, "coverage_by_gene.tsv")), delimiter="\t"))
    rows.sort(key=lambda r: float(r["raw_cov"]))
    genes = [r["gene"] for r in rows]; cov = [float(r["raw_cov"]) * 100 for r in rows]
    fig, ax = plt.subplots(figsize=(8, max(4, len(genes) * 0.16)))
    colors = ["#d62728" if c < 80 else "#ff7f0e" if c < 95 else "#2ca02c" for c in cov]
    ax.barh(genes, cov, color=colors)
    ax.axvline(95, ls="--", c="grey", lw=1)
    ax.set_xlabel("% of raw target covered"); ax.set_xlim(0, 100)
    ax.set_title("Per-gene coverage (raw); red<80, orange<95, green≥95")
    ax.tick_params(axis="y", labelsize=6)
    save(fig, figdir, "coverage_by_gene")


def fig_gc_tm(W, figdir):
    rows = list(csv.DictReader(open(os.path.join(W, "composition.tsv")), delimiter="\t"))
    gc_all = [float(r["gc"]) * 100 for r in rows]
    gc_pass = [float(r["gc"]) * 100 for r in rows if r["verdict"] == "PASS"]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(gc_all, bins=40, alpha=0.5, label=f"all spec-pass (n={len(gc_all)})", color="#7f7f7f")
    ax.hist(gc_pass, bins=40, alpha=0.7, label=f"composition-pass (n={len(gc_pass)})", color="#1f77b4")
    ax.axvline(25, ls="--", c="red", lw=1); ax.axvline(70, ls="--", c="red", lw=1)
    ax.set_xlabel("GC %"); ax.set_ylabel("baits"); ax.set_title("GC distribution (candidate vs passing)")
    ax.legend()
    save(fig, figdir, "gc_hist")

    tm = [float(r["tm"]) for r in rows if r["verdict"] == "PASS" and r["tm"] not in ("", "nan")]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(tm, bins=40, color="#ff7f0e")
    ax.axvline(60, ls="--", c="red", lw=1); ax.axvline(90, ls="--", c="red", lw=1)
    ax.set_xlabel("Tm (°C)"); ax.set_ylabel("baits"); ax.set_title("Tm distribution (passing baits)")
    save(fig, figdir, "tm_hist")


def fig_offtarget(W, figdir):
    s4 = load_stage(W, 4)
    dist = s4.get("offtarget_count_distribution", {})
    if not dist:
        return
    k = [int(x) for x in dist]; v = [dist[str(x)] for x in k]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar([str(x) for x in k], v, color="#d62728")
    ax.set_xlabel("off-target hit count"); ax.set_ylabel("baits (log)"); ax.set_yscale("log")
    ax.set_title("Off-target hit-count distribution (specificity)")
    save(fig, figdir, "offtarget_hist")


def fig_copynumber(W, figdir):
    s7 = load_stage(W, 7)
    dist = s7.get("copy_number_distribution", {})
    if not dist:
        return
    k = sorted(int(x) for x in dist)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar([str(x) for x in k], [dist[str(x)] for x in k], color="#9467bd")
    ax.set_xlabel("copy number"); ax.set_ylabel("unique baits")
    ax.set_title("Bait copy-number distribution (after boost/dedup)")
    save(fig, figdir, "copynumber_hist")


def fig_convergence(W, figdir):
    s6 = load_stage(W, 6)
    passes = s6.get("passes", [])
    if not passes:
        return
    xs = [p["pass"] for p in passes]
    ov = [p["coverage"].get("OVERALL", 0) * 100 for p in passes]
    ov.append(s6.get("final_coverage", {}).get("OVERALL", ov[-1] / 100) * 100)
    xs.append(xs[-1] + 1 if xs else 1)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(xs, ov, "-o", color="#2ca02c")
    ax.axhline(s6.get("coverage_target", 0.95) * 100, ls="--", c="grey", label="target")
    ax.set_xlabel("re-tiling pass"); ax.set_ylabel("mappable coverage %")
    ax.set_title(f"Coverage convergence (converged={s6.get('converged')})"); ax.legend()
    save(fig, figdir, "convergence")


def fig_per_class_summary(W, figdir):
    m = load_json(os.path.join(W, "metrics.json"))
    bc = m.get("coverage_by_class", {})
    final = defaultdict(int)
    for l in open(os.path.join(W, "baits.dedup.bed")):
        p = l.split("\t")[3].split("|")
        if len(p) >= 3:
            final[p[2]] += 1
    cls = [c for c in ["coding", "promoter", "regulatory", "cnv", "str", "anchor"] if c in bc or c in final]
    counts = [final.get(c, 0) for c in cls]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(cls, counts, color=[CLASS_COLORS.get(c, "#333") for c in cls])
    for i, c in enumerate(counts):
        ax.text(i, c, str(c), ha="center", va="bottom", fontsize=8)
    ax.set_ylabel("final unique baits"); ax.set_title("Final baits per design class")
    save(fig, figdir, "per_class_summary")


def fig_cost_tier(W, figdir):
    m = load_json(os.path.join(W, "metrics.json"))
    uniq = m.get("n_unique_baits") or 0
    tiers = m.get("pool_tiers", [12000, 92000, 244000])
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(["unique baits"], [uniq], color="#1f77b4", width=0.5)
    for t in tiers:
        ax.axhline(t, ls="--", c="grey", lw=1)
        ax.text(0.3, t, f"{t:,} tier", va="bottom", fontsize=8, color="grey")
    ax.text(0, uniq, f"{uniq:,}", ha="center", va="bottom", fontweight="bold")
    fit = m.get("pool_tier_fit")
    ax.set_ylabel("oligo count")
    ax.set_title(f"Unique baits vs pool tiers (fits {fit:,}; ${m.get('synthesis_pool_cost_usd')})"
                 if fit else "Unique baits vs pool tiers")
    save(fig, figdir, "cost_pool_tier")


def fig_coverage_tracks(W, figdir, per_page=6):
    raw = load_bed_by_gene(os.path.join(W, "targets.raw.bed"), 0)
    masked = load_bed_by_gene(os.path.join(W, "targets.masked.bed"), 0)
    baits = load_bed_by_gene(os.path.join(W, "baits.dedup.bed"), 1)   # CF_x|gene|class
    genes = sorted(raw.keys())
    pages = math.ceil(len(genes) / per_page)
    for pg in range(pages):
        sub = genes[pg * per_page:(pg + 1) * per_page]
        fig, axes = plt.subplots(len(sub), 1, figsize=(10, 1.7 * len(sub)), squeeze=False)
        for ax, gene in zip(axes[:, 0], sub):
            ivs = raw.get(gene, [])
            lo = min(s for _, s, e, _ in ivs); hi = max(e for _, s, e, _ in ivs)
            span = hi - lo or 1; pad = span * 0.02
            # raw targets by class (lane y=2)
            for ch, s, e, cls in ivs:
                ax.broken_barh([(s, e - s)], (1.6, 0.8), color=CLASS_COLORS.get(cls, "#333"), alpha=0.85)
            # masked-tileable (lane y=1) light
            for ch, s, e, cls in masked.get(gene, []):
                ax.broken_barh([(s, e - s)], (0.9, 0.5), color="#bbbbbb")
            # final baits (lane y=0)
            for ch, s, e, cls in baits.get(gene, []):
                ax.broken_barh([(s, e - s)], (0.2, 0.5), color="black")
            ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(0, 2.6)
            ax.set_yticks([0.45, 1.15, 2.0]); ax.set_yticklabels(["baits", "tileable", "target"], fontsize=7)
            ax.set_title(f"{gene}  ({ivs[0][0]}:{lo:,}-{hi:,})", fontsize=8, loc="left")
            ax.grid(False)
        present = sorted({c for g in sub for *_, c in raw.get(g, [])})
        fig.legend(handles=[Patch(color=CLASS_COLORS.get(c, "#333"), label=c) for c in present],
                   loc="lower center", ncol=len(present), fontsize=7, bbox_to_anchor=(0.5, -0.02))
        fig.suptitle(f"Per-gene coverage tracks (page {pg+1}/{pages})", y=1.0)
        save(fig, figdir, f"coverage_tracks_page{pg+1}")
    return pages


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--tracks-per-page", type=int, default=6)
    args = ap.parse_args()
    os.makedirs(args.figdir, exist_ok=True)
    W = args.workdir
    fig_coverage_by_class(W, args.figdir)
    fig_coverage_by_gene(W, args.figdir)
    fig_gc_tm(W, args.figdir)
    fig_offtarget(W, args.figdir)
    fig_copynumber(W, args.figdir)
    fig_convergence(W, args.figdir)
    fig_per_class_summary(W, args.figdir)
    fig_cost_tier(W, args.figdir)
    pages = fig_coverage_tracks(W, args.figdir, args.tracks_per_page)
    import sys
    sys.stderr.write(f"[make_figures] wrote figures (PNG+SVG) to {args.figdir} "
                     f"(incl. {pages} coverage-track pages)\n")


if __name__ == "__main__":
    main()
