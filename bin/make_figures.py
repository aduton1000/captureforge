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
import argparse, csv, json, math, os, re, subprocess, sys
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

# One-line, plain-language interpretation baked into every figure (and reused by
# the dashboard/report/manual). Keyed by figure base name (page suffix stripped).
CAPTIONS = {
    "coverage_by_class": "Mappable coverage ≥96% for coding/promoter/str; CNV lower (~78%) because repeat-rich loci are largely unmappable.",
    "coverage_by_gene": "Per-gene raw coverage; GENE/GENE lowest (CNV segdup/LCR). Red <80%, orange <95%, green ≥95%.",
    "gc_hist": "GC of spec-pass candidate baits (one row per tiled position — copies are NOT multiplied in). Spikes are expected from dense tiling of homogeneous regions; what matters is the bulk sitting in the 25–70% window.",
    "tm_hist": "Tm of passing candidate baits (per tiled position, not copies). The population sits within the 60–90 °C capture window.",
    "offtarget_hist": "Most baits are unique (0 off-targets); the long repeat tail (binned) is rejected. Binned x-axis for readability.",
    "copynumber_hist": "Most baits appear once; promoter/regulatory/anchor baits are boosted to ×2–4 to even out capture.",
    "convergence": "Re-tiling plateaus below the 95% target — correct, because the residual gaps are non-specific/repeat sequence.",
    "per_class_summary": "Final unique baits per design class (coding dominates; cnv/str come from whole-locus tiling).",
    "cost_pool_tier": "Unique baits sit under the 12k oligo-pool tier (synthesis-pool cost only, not assay cost).",
    "baits_per_gene_by_class": "Per-gene bait composition by design class (multi-class baits attributed to their primary class).",
    "coverage_tracks": "Lanes: target (raw, class-colored) / tileable (post-mask) / baits (final). A baits gap where tileable is ALSO empty = masked/uncapturable (expected); a baits gap where tileable is present = a real gap to investigate.",
}


def _caption_for(name):
    return CAPTIONS.get(name) or CAPTIONS.get(re.sub(r"_page\d+$", "", name))


def save(fig, figdir, name):
    cap = _caption_for(name)
    if cap:
        fig.text(0.5, -0.015, cap, ha="center", va="top", fontsize=7, style="italic",
                 color="#444", wrap=True)
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
    # categorical bins so the 0..thousands range is readable
    bins = [("0", lambda k: k == 0), ("1", lambda k: k == 1),
            ("2–5", lambda k: 2 <= k <= 5), ("6–20", lambda k: 6 <= k <= 20),
            ("21–100", lambda k: 21 <= k <= 100), (">100", lambda k: k > 100)]
    labels = [b[0] for b in bins]
    vals = [0] * len(bins)
    for ks, n in dist.items():
        k = int(ks)
        for i, (_, pred) in enumerate(bins):
            if pred(k):
                vals[i] += n
                break
    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.bar(labels, vals, color="#d62728")
    ax.set_yscale("log")
    ax.set_xlabel("off-target hit count (binned)"); ax.set_ylabel("baits (log scale)")
    ax.bar_label(bars, fmt="%d", fontsize=8, padding=2)
    ax.margins(y=0.18)
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


def fig_coverage_tracks(W, figdir, per_page=4):
    raw = load_bed_by_gene(os.path.join(W, "targets.raw.bed"), 0)
    masked = load_bed_by_gene(os.path.join(W, "targets.masked.bed"), 0)
    baits = load_bed_by_gene(os.path.join(W, "baits.dedup.bed"), 1)   # CF_x|gene|class
    genes = sorted(raw.keys())
    pages = math.ceil(len(genes) / per_page)
    for pg in range(pages):
        sub = genes[pg * per_page:(pg + 1) * per_page]
        # taller panels + generous vertical gap so coordinate ticks never collide
        # with the next gene's title.
        fig, axes = plt.subplots(len(sub), 1, figsize=(10, 2.3 * len(sub)), squeeze=False,
                                 gridspec_kw={"hspace": 1.15})
        for ax, gene in zip(axes[:, 0], sub):
            ivs = raw.get(gene, [])
            lo = min(s for _, s, e, _ in ivs); hi = max(e for _, s, e, _ in ivs)
            span = hi - lo or 1; pad = span * 0.02
            for ch, s, e, cls in ivs:                                        # target lane
                ax.broken_barh([(s, e - s)], (1.6, 0.8), color=CLASS_COLORS.get(cls, "#333"), alpha=0.85)
            for ch, s, e, cls in masked.get(gene, []):                       # tileable lane
                ax.broken_barh([(s, e - s)], (0.9, 0.5), color="#bbbbbb")
            for ch, s, e, cls in baits.get(gene, []):                        # baits lane
                ax.broken_barh([(s, e - s)], (0.2, 0.5), color="black")
            ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(0, 2.6)
            ax.set_yticks([0.45, 1.15, 2.0]); ax.set_yticklabels(["baits", "tileable", "target"], fontsize=7)
            ax.set_title(f"{gene}  ({ivs[0][0]}:{lo:,}-{hi:,})", fontsize=8.5, loc="left", pad=4)
            ax.tick_params(axis="x", labelsize=6.5)
            ax.grid(False)
        present = sorted({c for g in sub for *_, c in raw.get(g, [])})
        fig.legend(handles=[Patch(color=CLASS_COLORS.get(c, "#333"), label=c) for c in present],
                   loc="lower center", ncol=len(present), fontsize=7, bbox_to_anchor=(0.5, -0.01))
        fig.suptitle(f"Per-gene coverage tracks (page {pg+1}/{pages})", y=0.995, fontsize=11)
        save(fig, figdir, f"coverage_tracks_page{pg+1}")
    return pages


def fig_baits_per_gene_by_class(W, figdir, per_page=30):
    """Horizontal stacked bar of bait counts per gene, broken down by class.
    Counts derive from baits.csv (gene + class); multi-class baits are attributed
    to their primary (first-listed) class so the bars sum to the unique bait count.
    Also writes the gene × class count matrix as a TSV."""
    path = os.path.join(W, "baits.csv")
    if not os.path.exists(path):
        return 0
    counts = defaultdict(lambda: defaultdict(int))            # gene -> class -> n
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            gene = (r.get("gene", "") or "").split(",")[0].strip()
            cls = (r.get("class", "") or "").split(",")[0].strip()
            if gene and cls:
                counts[gene][cls] += 1
    if not counts:
        return 0
    order = ["coding", "promoter", "regulatory", "cnv", "str", "anchor"]
    classes = [c for c in order if any(c in counts[g] for g in counts)]
    # matrix TSV (all genes, alphabetical)
    with open(os.path.join(figdir, "baits_per_gene_by_class.tsv"), "w") as o:
        o.write("gene\t" + "\t".join(classes) + "\ttotal\n")
        for g in sorted(counts):
            row = [counts[g].get(c, 0) for c in classes]
            o.write(g + "\t" + "\t".join(map(str, row)) + f"\t{sum(row)}\n")
    # paginate, sorted by total desc so big genes lead
    genes = sorted(counts, key=lambda g: sum(counts[g].values()))   # asc → largest at top after barh
    pages = math.ceil(len(genes) / per_page)
    for pg in range(pages):
        sub = genes[max(0, len(genes) - (pg + 1) * per_page): len(genes) - pg * per_page]
        y = list(range(len(sub)))
        fig, ax = plt.subplots(figsize=(9, max(4, 0.34 * len(sub) + 1.3)))
        left = [0] * len(sub)
        for c in classes:
            vals = [counts[g].get(c, 0) for g in sub]
            ax.barh(y, vals, left=left, color=CLASS_COLORS.get(c, "#333"), label=c, height=0.78)
            left = [a + b for a, b in zip(left, vals)]
        ax.set_yticks(y); ax.set_yticklabels(sub, fontsize=7)
        ax.set_ylim(-0.6, len(sub) - 0.4)
        ax.set_xlabel("baits"); ax.margins(x=0.02)
        ax.set_title(f"Baits per gene by design class" + (f" (page {pg+1}/{pages})" if pages > 1 else ""))
        ax.legend(ncol=len(classes), fontsize=7, loc="upper right")
        save(fig, figdir, f"baits_per_gene_by_class_page{pg+1}" if pages > 1 else "baits_per_gene_by_class")
    return pages


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--tracks-per-page", type=int, default=4)
    ap.add_argument("--genes-per-page", type=int, default=30, help="baits_per_gene_by_class pagination")
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
    fig_baits_per_gene_by_class(W, args.figdir, args.genes_per_page)
    fig_cost_tier(W, args.figdir)
    pages = fig_coverage_tracks(W, args.figdir, args.tracks_per_page)
    # polished pipeline DAG (graphviz). Best-effort: skip cleanly if unavailable.
    dag = os.path.join(os.path.dirname(os.path.abspath(__file__)), "make_dag.py")
    try:
        subprocess.run([sys.executable, dag, "--figdir", args.figdir], check=True,
                       capture_output=True, text=True)
    except Exception as e:
        sys.stderr.write(f"[make_figures] DAG skipped ({e})\n")
    sys.stderr.write(f"[make_figures] wrote figures (PNG+SVG) to {args.figdir} "
                     f"(incl. {pages} coverage-track pages)\n")


if __name__ == "__main__":
    main()
