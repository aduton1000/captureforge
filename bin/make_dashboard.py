#!/usr/bin/env python3
"""
CaptureForge — visualization : make_dashboard.py

Build a single self-contained HTML QC dashboard: embeds every figure (base64
PNG) and the key metrics tables from metrics.json. No external assets.
"""
from __future__ import annotations
import argparse, base64, glob, json, os, re, sys
from pathlib import Path
from jinja2 import Template

# Plain-language one-liner per figure (kept in sync with bin/make_figures.py).
CAPTIONS = {
    "pipeline_dag": "The 10-stage CaptureForge flow, colored by phase: target definition → selection → assembly → QC & export.",
    "coverage_by_class": "Mappable coverage ≥96% for coding/promoter/str; CNV lower (~78%) because repeat-rich loci are largely unmappable.",
    "coverage_by_gene": "Per-gene raw coverage; GENE/GENE lowest (CNV segdup/LCR). Red <80%, orange <95%, green ≥95%.",
    "baits_per_gene_by_class": "Per-gene bait composition by design class (multi-class baits attributed to their primary class).",
    "gc_hist": "GC of spec-pass candidate baits — one row per tiled position, copies NOT multiplied in. Spikes are expected from dense tiling; what matters is the bulk inside the 25–70% window.",
    "tm_hist": "Tm of passing candidate baits (per position, not copies); the population sits within the 60–90 °C capture window.",
    "offtarget_hist": "Most baits are unique (0 off-targets); the long repeat tail (binned) is rejected.",
    "copynumber_hist": "Most baits appear once; promoter/regulatory/anchor baits are boosted ×2–4 to even out capture.",
    "per_class_summary": "Final unique baits per design class (coding dominates; cnv/str from whole-locus tiling).",
    "convergence": "Re-tiling plateaus below the 95% target — correct, because residual gaps are non-specific/repeat.",
    "cost_pool_tier": "Unique baits sit under the 12k oligo-pool tier (synthesis-pool cost only).",
    "bait_construct": "What the synthesized oligo looks like for each format (segments scaled by nt). ivt_rna = 5′ handle + T7 + bait + 3′ handle; dna_pool / ready_to_use = bait core only. The configured format is marked.",
    "coverage_tracks": "Lanes: target (raw) / exons (coding exons ±15 bp) / tileable (post-mask) / baits (final). A baits gap where tileable is also empty = masked/uncapturable (expected); a gap where tileable is present = a real gap. The exon lane shows whether any exon is missed.",
}

def caption_for(name):
    return CAPTIONS.get(name) or CAPTIONS.get(re.sub(r"_page\d+$", "", name)) or ""

TEMPLATE = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>CaptureForge QC — {{ panel }}</title>
<style>
 body{font-family:-apple-system,Segoe UI,Roboto,Arial,sans-serif;margin:0;background:#f6f7f9;color:#1b1f24}
 header{background:#0d3b66;color:#fff;padding:18px 28px}
 header h1{margin:0;font-size:21px} header .sub{opacity:.85;font-size:13px;margin-top:4px}
 main{max-width:1100px;margin:0 auto;padding:22px}
 .cards{display:flex;flex-wrap:wrap;gap:12px;margin:8px 0 22px}
 .card{background:#fff;border:1px solid #e2e6ea;border-radius:10px;padding:14px 18px;min-width:150px;flex:1}
 .card .v{font-size:24px;font-weight:700;color:#0d3b66} .card .k{font-size:12px;color:#5b6470;margin-top:2px}
 .card.warn .v{color:#b3261e}
 section{background:#fff;border:1px solid #e2e6ea;border-radius:10px;padding:16px 20px;margin-bottom:18px}
 h2{font-size:16px;border-bottom:2px solid #0d3b66;padding-bottom:6px;margin-top:0}
 table{border-collapse:collapse;width:100%;font-size:13px} th,td{border:1px solid #e2e6ea;padding:5px 9px;text-align:right}
 th:first-child,td:first-child{text-align:left} th{background:#eef2f6}
 img{max-width:100%;border:1px solid #e8ebee;border-radius:6px;margin:8px 0}
 .fig{margin-bottom:14px} .fig .cap{font-size:12px;color:#5b6470;margin-top:3px}
 details{margin:6px 0} summary{cursor:pointer;font-weight:600}
 .flag{background:#fdecea;border-left:4px solid #b3261e;padding:8px 12px;border-radius:5px;margin:6px 0;font-size:13px}
 .ok{background:#e7f4ea;border-left:4px solid #1e7d32;padding:8px 12px;border-radius:5px;margin:6px 0;font-size:13px}
</style></head><body>
<header><h1>CaptureForge — Hybridization-Capture Probe Design QC</h1>
<div class="sub">panel: <b>{{ panel }}</b> &nbsp;•&nbsp; generated {{ generated }} &nbsp;•&nbsp; bait_format: {{ bait_format }}</div></header>
<main>
 <div class="cards">
  <div class="card"><div class="v">{{ m.n_unique_baits or "—" }}</div><div class="k">unique baits</div></div>
  <div class="card"><div class="v">{{ (m.coverage_overall.mappable_cov*100)|round(1) }}%</div><div class="k">mappable coverage</div></div>
  <div class="card"><div class="v">{{ (m.coverage_overall.raw_cov*100)|round(1) }}%</div><div class="k">raw coverage</div></div>
  <div class="card"><div class="v">${{ m.synthesis_pool_cost_usd }}</div><div class="k">pool synthesis cost</div></div>
  <div class="card {{ 'warn' if m.tripwire_exceeds_tier else '' }}"><div class="v">{{ '{:,}'.format(m.pool_tier_fit) if m.pool_tier_fit else 'OVER' }}</div><div class="k">pool tier fit</div></div>
 </div>

 {% if m.tripwire_exceeds_tier %}<div class="flag">⚠ Tripwire: unique baits exceed the smallest pool tier — decide trim-vs-pay.</div>
 {% else %}<div class="ok">✓ Fits the {{ '{:,}'.format(m.pool_tier_fit) }} oligo-pool tier (pool synthesis cost only).</div>{% endif %}
 <div class="ok">Coverage is reported as <b>mappable</b> (over tileable sequence) and <b>raw</b> (over nominal target). The
  raw↔mappable gap is sequence intentionally removed by masking (repeats/segdups/STR cores/Alu) — uncapturable by design.</div>

 {% if has_dag %}<section><h2>Pipeline</h2>
  {{ figure('pipeline_dag') }}
 </section>{% endif %}

 <section><h2>Coverage</h2>
  {% for f in ['coverage_by_class','coverage_by_gene'] %}{{ figure(f) }}{% endfor %}
  <table><tr><th>class</th><th>raw cov</th><th>mappable cov</th><th>covered bp</th><th>raw bp</th></tr>
  {% for c,v in m.coverage_by_class.items() %}<tr><td>{{c}}</td><td>{{(v.raw_cov*100)|round(1)}}%</td>
   <td>{{(v.mappable_cov*100)|round(1)}}%</td><td>{{v.covered_bp}}</td><td>{{v.raw_bp}}</td></tr>{% endfor %}</table>
 </section>

 <section><h2>CNV read-depth callability</h2>
  <table><tr><th>gene</th><th>bins</th><th>median gap</th><th>max gap</th><th>CV</th><th>depth-callable</th></tr>
  {% for g,v in m.cnv_bin_spacing.items() %}<tr><td>{{g}}</td><td>{{v.n_bins}}</td><td>{{v.median_gap}}</td>
   <td>{{v.max_gap}}</td><td>{{v.cv_gap}}</td><td>{{ '✓' if v.depth_callable else '✗ (gaps at segdups)' }}</td></tr>{% endfor %}</table>
  <div class="cap">Breakpoint-targeted CNV design is a future enhancement (needs curated breakpoints).</div>
 </section>

 <section><h2>Bait composition & specificity</h2>
  {% for f in ['gc_hist','tm_hist','offtarget_hist','copynumber_hist'] %}{{ figure(f) }}{% endfor %}
 </section>

 <section><h2>Per-gene bait composition</h2>
  {% for f in pergene_pages %}{{ figure(f) }}{% endfor %}
 </section>

 <section><h2>Tiling / convergence / pool</h2>
  {% for f in ['per_class_summary','convergence','cost_pool_tier'] %}{{ figure(f) }}{% endfor %}
 </section>

 <section><h2>Bait construct (synthesis)</h2>
  {{ figure('bait_construct') }}
 </section>

 <section><h2>Per-gene coverage tracks</h2>
  {% for f in track_pages %}{{ figure(f) }}{% endfor %}
 </section>

 <section><h2>Raw metrics</h2>
  <details><summary>metrics.json</summary><pre style="font-size:11px;overflow:auto">{{ metrics_json }}</pre></details>
 </section>
</main></body></html>"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metrics", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--panel", default="captureforge_panel")
    ap.add_argument("--bait-format", default="ivt_rna")
    ap.add_argument("--generated", default="")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    m = json.loads(Path(args.metrics).read_text())

    def embed(name):
        p = os.path.join(args.figdir, f"{name}.png")
        if not os.path.exists(p):
            return None
        return base64.b64encode(Path(p).read_bytes()).decode()

    def figure(name):
        b64 = embed(name)
        if not b64:
            return ""
        cap = caption_for(name)
        cap_html = f'<div class="cap"><b>{cap}</b></div>' if cap else ""
        return (f'<div class="fig"><img src="data:image/png;base64,{b64}">{cap_html}'
                f'<div class="cap" style="opacity:.6">{name}.png (SVG also in figures/)</div></div>')

    def pages(prefix):
        return sorted(os.path.basename(p)[:-4]
                      for p in glob.glob(os.path.join(args.figdir, f"{prefix}*.png")))

    track_pages = pages("coverage_tracks_page")
    pergene_pages = (["baits_per_gene_by_class"]
                     if os.path.exists(os.path.join(args.figdir, "baits_per_gene_by_class.png"))
                     else pages("baits_per_gene_by_class_page"))
    has_dag = os.path.exists(os.path.join(args.figdir, "pipeline_dag.png"))

    html = Template(TEMPLATE).render(
        m=m, panel=args.panel, bait_format=args.bait_format,
        generated=args.generated, metrics_json=json.dumps(m, indent=2),
        figure=figure, track_pages=track_pages, pergene_pages=pergene_pages, has_dag=has_dag)
    Path(args.output).write_text(html)
    sys.stderr.write(f"[make_dashboard] wrote self-contained dashboard → {args.output} "
                     f"({len(html)//1024} KB, {len(track_pages)} track pages)\n")


if __name__ == "__main__":
    main()
