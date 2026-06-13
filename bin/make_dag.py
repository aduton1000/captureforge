#!/usr/bin/env python3
"""
CaptureForge — visualization : make_dag.py

Render a polished, color-coded 10-stage pipeline DAG with graphviz: nodes grouped
by phase, each labelled with the stage name + its primary tool, clean data-flow
edges, and a legend. Writes pipeline_dag.png (high-res) + pipeline_dag.svg.

Falls back gracefully (non-zero exit handled by caller) if graphviz is absent.
"""
from __future__ import annotations
import argparse, os, sys

# phase palette (distinct but harmonious)
PHASE = {
    "io":       ("#cfd8dc", "#455a64"),   # fill, border
    "define":   ("#bbdefb", "#1565c0"),
    "filter":   ("#ffe0b2", "#e65100"),
    "assemble": ("#c8e6c9", "#2e7d32"),
    "report":   ("#e1bee7", "#6a1b9a"),
}
PHASE_LABEL = {
    "define": "Target definition (Stages 1–3)",
    "filter": "Selection (Stages 4–6)",
    "assemble": "Assembly (Stages 7–8)",
    "report": "QC & export (Stages 9–10)",
}

# id, label, phase
NODES = [
    ("in",  "Inputs\\ngene list · GTF · no-alt genome", "io"),
    ("s1",  "Stage 1 · Target definition\\ndefine_targets.py (GTF / MANE)", "define"),
    ("s2",  "Stage 2 · Masking\\napply_masking.py (dustmasker + rmsk/simpleRepeat)", "define"),
    ("s3",  "Stage 3 · Tiling\\ntile_candidates.py (bedtools)", "define"),
    ("s4",  "Stage 4 · Specificity\\nspecificity_filter.py (blastn / bowtie2)", "filter"),
    ("s5",  "Stage 5 · Composition\\ncomposition_filter.py (primer3)", "filter"),
    ("s6",  "Stage 6 · Re-tiling\\ngap_retile.py", "filter"),
    ("s7",  "Stage 7 · Boost / dedup\\nboost_dedup.py", "assemble"),
    ("s8",  "Stage 8 · Synthesis scaffolding\\nscaffold_synthesis.py", "assemble"),
    ("s9",  "Stage 9 · QC aggregation\\naggregate_qc.py", "report"),
    ("s10", "Stage 10 · Export\\nexport_package.py", "report"),
    ("out", "Outputs\\nbaits.fasta/csv · final_covered_targets.bed · qc_dashboard.html", "io"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--name", default="pipeline_dag")
    args = ap.parse_args()
    os.makedirs(args.figdir, exist_ok=True)

    try:
        from graphviz import Digraph
    except ImportError:
        sys.stderr.write("[make_dag] python-graphviz not installed; skipping DAG\n")
        sys.exit(2)

    g = Digraph("captureforge", engine="dot")
    g.attr(rankdir="TB", bgcolor="white", fontname="Helvetica",
           nodesep="0.35", ranksep="0.5", dpi="200", label=" ", labelloc="t")
    g.attr("node", shape="box", style="rounded,filled", fontname="Helvetica",
           fontsize="11", margin="0.20,0.12", penwidth="1.4")
    g.attr("edge", color="#90a4ae", arrowsize="0.8", penwidth="1.6")

    for nid, label, phase in NODES:
        fill, border = PHASE[phase]
        g.node(nid, label, fillcolor=fill, color=border)

    chain = [n[0] for n in NODES]
    for a, b in zip(chain, chain[1:]):
        g.edge(a, b)

    # legend cluster (phase swatches), kept off to the side
    with g.subgraph(name="cluster_legend") as lg:
        lg.attr(label="Phase", fontname="Helvetica", fontsize="11",
                color="#cfd8dc", style="rounded", margin="8")
        lg.attr("node", shape="box", style="filled", fontsize="10", width="2.6",
                height="0.3", fixedsize="false", penwidth="1.2")
        prev = None
        for ph in ("define", "filter", "assemble", "report"):
            fill, border = PHASE[ph]
            nid = f"lg_{ph}"
            lg.node(nid, PHASE_LABEL[ph], fillcolor=fill, color=border)
            if prev:
                lg.edge(prev, nid, style="invis")
            prev = nid

    base = os.path.join(args.figdir, args.name)
    g.format = "png"; g.render(base, cleanup=True)
    g.format = "svg"; g.render(base, cleanup=True)
    sys.stderr.write(f"[make_dag] wrote {base}.png and {base}.svg\n")


if __name__ == "__main__":
    main()
