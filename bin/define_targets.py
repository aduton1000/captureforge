#!/usr/bin/env python3
"""
CaptureForge — Stage 1 : define_targets.py

Turn the normalized target table (parse_gene_list.py output) into a labelled,
design-classed BED of capture targets, using a local Ensembl GTF for canonical
transcript structure. Every interval is tagged `gene|class|subtype` so the
masking and tiling stages can apply per-class policy.

Per design class (see assets/gene_list_schema.md and docs/design_decisions.md):
  coding      : canonical-transcript exons (incl. UTRs) ± splice_flank
  promoter    : strand-aware TSS window  (TSS-up .. TSS+down)
  regulatory  : Ensembl Regulatory / cCRE features overlapping the locus (if provided)
  cnv         : whole gene locus (tiled EVENLY downstream)
  str         : whole gene locus (repeat core removed at masking → unique flanks only)
  anchor      : rsID ± anchor_flank (needs dbSNP VCF or --rsid-bed; else skipped+flagged)

Usage:
  define_targets.py --targets targets.tsv --gtf ann.gtf.gz --fai genome.fa.fai \
     [--canonical mane] [--splice-flank 15] [--promoter-up 2000 --promoter-down 200] \
     [--anchor-flank 150] [--dbsnp dbsnp.vcf.gz] [--regulatory reg.gff.gz] \
     -o targets.raw.bed --summary stage1_summary.json
"""
from __future__ import annotations
import argparse, csv, gzip, json, re, sys
from collections import defaultdict
from pathlib import Path

ATTR = {k: re.compile(rf'{k} "([^"]+)"') for k in ("gene_id", "transcript_id", "gene_name")}
TAG_RE = re.compile(r'tag "([^"]+)"')


def openmaybe(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def load_fai(fai):
    sizes = {}
    with open(fai) as fh:
        for line in fh:
            f = line.split("\t")
            sizes[f[0]] = int(f[1])
    return sizes


def parse_gtf_for_genes(gtf, wanted_ids):
    """Collect transcript + exon structure for the wanted gene_ids only."""
    tx = {}            # tid -> dict(gene_id, chrom, start, end, strand, tags:set, exons:list)
    with openmaybe(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            mgid = ATTR["gene_id"].search(f[8])
            if not mgid or mgid.group(1) not in wanted_ids:
                continue
            feat = f[2]
            if feat == "transcript":
                mtid = ATTR["transcript_id"].search(f[8])
                if not mtid:
                    continue
                tx[mtid.group(1)] = {
                    "gene_id": mgid.group(1), "chrom": f[0],
                    "start": int(f[3]), "end": int(f[4]), "strand": f[6],
                    "tags": set(TAG_RE.findall(f[8])), "exons": [],
                }
            elif feat == "exon":
                mtid = ATTR["transcript_id"].search(f[8])
                if mtid and mtid.group(1) in tx:
                    tx[mtid.group(1)]["exons"].append((int(f[3]), int(f[4])))
    return tx


def pick_transcript(transcripts, mode):
    """transcripts: list of tx dicts for one gene. Return chosen tx dict."""
    if not transcripts:
        return None
    if mode == "mane":
        for t in transcripts:
            if "MANE_Select" in t["tags"]:
                return t
    if mode in ("mane", "ensembl_canonical"):
        for t in transcripts:
            if "Ensembl_canonical" in t["tags"]:
                return t
    # longest by summed exon length
    return max(transcripts, key=lambda t: sum(e - s + 1 for s, e in t["exons"]) or (t["end"] - t["start"]))


def clip(start, end, chrom, sizes):
    """0-based half-open clip to chromosome bounds."""
    lo = max(0, start)
    hi = min(sizes.get(chrom, end), end)
    return (lo, hi) if hi > lo else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--fai", required=True, help="genome .fai for chrom-size clipping")
    ap.add_argument("--canonical", default="mane", choices=["mane", "ensembl_canonical", "longest"])
    ap.add_argument("--splice-flank", type=int, default=15)
    ap.add_argument("--promoter-up", type=int, default=2000)
    ap.add_argument("--promoter-down", type=int, default=200)
    ap.add_argument("--anchor-flank", type=int, default=150)
    ap.add_argument("--dbsnp", default=None, help="dbSNP VCF(.gz, tabix-indexed) for rsID anchors")
    ap.add_argument("--rsid-bed", default=None, help="optional rsID->BED (chrom start end rsid) fallback")
    ap.add_argument("--regulatory", default=None, help="Ensembl Regulatory / cCRE GFF(.gz)")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--summary", default=None)
    args = ap.parse_args()

    sizes = load_fai(args.fai)

    # read normalized target table
    rows = list(csv.DictReader(open(args.targets), delimiter="\t"))
    if not rows:
        sys.exit("ERROR: empty targets table")
    wanted_ids = {r["ensembl_id"] for r in rows}
    sys.stderr.write(f"[define_targets] genes={len(rows)}  parsing GTF…\n")
    tx_all = parse_gtf_for_genes(args.gtf, wanted_ids)
    by_gene = defaultdict(list)
    for t in tx_all.values():
        by_gene[t["gene_id"]].append(t)

    # optional rsID coordinate lookup
    rsid_coord = {}
    if args.rsid_bed:
        for line in open(args.rsid_bed):
            f = line.split()
            if len(f) >= 4:
                rsid_coord[f[3].lower()] = (f[0], int(f[1]), int(f[2]))
    dbsnp = None
    if args.dbsnp:
        try:
            import pysam
            dbsnp = pysam.TabixFile(args.dbsnp)
        except Exception as e:
            sys.stderr.write(f"[define_targets] WARN cannot open dbSNP ({e}); anchors via VCF disabled\n")

    out = open(args.output, "w")
    n_intervals = 0
    bp_by_gene_class = defaultdict(int)
    n_by_class = defaultdict(int)
    anchors_missing = []
    regulatory_features = load_regulatory(args.regulatory) if args.regulatory else None
    genes_no_tx = []

    def emit(chrom, start, end, gene, cls, subtype, strand):
        nonlocal n_intervals
        c = clip(start, end, chrom, sizes)
        if not c:
            return
        s, e = c
        out.write(f"{chrom}\t{s}\t{e}\t{gene}|{cls}|{subtype}\t0\t{strand}\n")
        n_intervals += 1
        bp_by_gene_class[(gene, cls)] += e - s
        n_by_class[cls] += 1

    for r in rows:
        gene, gid = r["gene"], r["ensembl_id"]
        classes = r["design_classes"].split("|") if r["design_classes"] else ["coding"]
        anchors = [a for a in r.get("anchors", "").split(",") if a]
        extra = [x for x in r.get("extra_regions", "").split(";") if x]
        t = pick_transcript(by_gene.get(gid, []), args.canonical)
        if not t:
            genes_no_tx.append(gene)
            continue
        chrom, strand = t["chrom"], t["strand"]
        gstart, gend = int(r["gene_start"]), int(r["gene_end"])

        if "coding" in classes:
            for i, (es, ee) in enumerate(sorted(t["exons"]), 1):
                emit(chrom, es - 1 - args.splice_flank, ee + args.splice_flank,
                     gene, "coding", f"exon{i}", strand)
        if "promoter" in classes:
            tss = t["start"] if strand == "+" else t["end"]
            if strand == "+":
                ps, pe = tss - args.promoter_up, tss + args.promoter_down
            else:
                ps, pe = tss - args.promoter_down, tss + args.promoter_up
            emit(chrom, ps - 1, pe, gene, "promoter", "tss", strand)
        if "regulatory" in classes and regulatory_features is not None:
            for (rs, re_, rtype) in regulatory_features.get(chrom, []):
                if rs <= gend and re_ >= gstart:        # overlaps locus
                    emit(chrom, rs, re_, gene, "regulatory", rtype, strand)
        if "cnv" in classes:
            emit(chrom, gstart - 1, gend, gene, "cnv", "locus", strand)
        if "str" in classes:
            emit(chrom, gstart - 1, gend, gene, "str", "flank", strand)
        if "anchor" in classes:
            for rs in anchors:
                coord = rsid_coord.get(rs.lower())
                if not coord and dbsnp is not None:
                    coord = lookup_dbsnp(dbsnp, rs, chrom, gstart, gend)
                if coord:
                    c, pos0, _ = coord
                    emit(c, pos0 - args.anchor_flank, pos0 + 1 + args.anchor_flank,
                         gene, "anchor", rs, strand)
                else:
                    anchors_missing.append(f"{gene}:{rs}")
        for x in extra:
            m = re.match(r"([^:]+):(\d+)-(\d+)", x)
            if m:
                emit(m.group(1), int(m.group(2)) - 1, int(m.group(3)), gene, "extra", "region", strand)

    out.close()

    summary = {
        "n_genes": len(rows),
        "n_intervals": n_intervals,
        "intervals_by_class": dict(n_by_class),
        "bp_by_class": {c: sum(bp for (g, cc), bp in bp_by_gene_class.items() if cc == c)
                        for c in n_by_class},
        "bp_by_gene": {g: sum(bp for (gg, cc), bp in bp_by_gene_class.items() if gg == g)
                       for g in {g for (g, _) in bp_by_gene_class}},
        "genes_without_transcript": genes_no_tx,
        "anchors_unresolved": anchors_missing,
        "canonical_source": args.canonical,
        "regulatory_used": bool(regulatory_features),
        "dbsnp_used": bool(dbsnp),
    }
    # outlier-footprint genes (>3x median bp)
    bpg = summary["bp_by_gene"]
    if bpg:
        vals = sorted(bpg.values()); med = vals[len(vals) // 2]
        summary["outlier_footprint_genes"] = {g: v for g, v in bpg.items() if med and v > 3 * med}
    if args.summary:
        Path(args.summary).write_text(json.dumps(summary, indent=2))
    sys.stderr.write(f"[define_targets] wrote {n_intervals} intervals  by_class={dict(n_by_class)}\n")
    if genes_no_tx:
        sys.stderr.write(f"[define_targets] WARN no transcript for: {genes_no_tx}\n")
    if anchors_missing:
        sys.stderr.write(f"[define_targets] NOTE {len(anchors_missing)} rsID anchors unresolved "
                         f"(no dbSNP/rsid-bed): {anchors_missing}\n")


def load_regulatory(gff):
    feats = defaultdict(list)
    with openmaybe(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            feats[f[0]].append((int(f[3]) - 1, int(f[4]), (f[2] or "reg")[:20]))
    return feats


def lookup_dbsnp(tbx, rsid, chrom, gstart, gend):
    """Search dbSNP VCF near the gene for the rsID; return (chrom, pos0, pos0)."""
    want = rsid.lower()
    try:
        for rec in tbx.fetch(chrom, max(0, gstart - 200000), gend + 200000):
            f = rec.split("\t")
            if len(f) >= 3 and f[2].lower() == want:
                return (f[0], int(f[1]) - 1, int(f[1]) - 1)
    except Exception:
        return None
    return None


if __name__ == "__main__":
    main()
