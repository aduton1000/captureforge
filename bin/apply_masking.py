#!/usr/bin/env python3
"""
CaptureForge — Stage 2 : apply_masking.py

Apply per-design-class masking policy to the raw targets and report maskability.

Masks:
  low_complexity : computed with dustmasker over the target sequence (always)
  hard_mask      : RepeatMasker + segmental dups (optional --rmsk/--segdups)

Per-class policy:
  coding / cnv         : subtract hard_mask + low_complexity   (standard exclude)
  promoter / regulatory: subtract hard_mask only               (RELAX low-complexity)
  str                  : subtract hard_mask + low_complexity    → unique FLANKS only;
                         the removed low-complexity core is written to *.str_core.bed
                         and flagged for a separate assay
  anchor / extra       : minimal masking (kept as-is)

Outputs: targets.masked.bed, targets.lost_to_mask.bed, <prefix>.str_core.bed, summary JSON.
"""
from __future__ import annotations
import argparse, json, os, re, subprocess, sys, tempfile
from collections import defaultdict
from pathlib import Path


def run(cmd, **kw):
    return subprocess.run(cmd, check=True, text=True, capture_output=True, **kw)


def read_bed(path):
    """Return dict: class -> dict: chrom -> list[(start,end,name,strand)]."""
    by_class = defaultdict(lambda: defaultdict(list))
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            chrom, s, e, name = f[0], int(f[1]), int(f[2]), f[3]
            strand = f[5] if len(f) > 5 else "+"
            cls = name.split("|")[1] if "|" in name else "coding"
            by_class[cls][chrom].append((s, e, name, strand))
    return by_class


def merge_intervals(ivs):
    """ivs: list[(s,e)] -> merged sorted list."""
    if not ivs:
        return []
    ivs = sorted(ivs)
    out = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def subtract(intervals, masks):
    """Subtract masks (list[(s,e)]) from each (s,e,name,strand). Returns kept + removed."""
    masks = merge_intervals(masks)
    kept, removed = [], []
    for s, e, name, strand in intervals:
        cur = [(s, e)]
        for ms, me in masks:
            nxt = []
            for cs, ce in cur:
                if me <= cs or ms >= ce:          # no overlap
                    nxt.append((cs, ce)); continue
                if ms > cs:
                    nxt.append((cs, ms))
                if me < ce:
                    nxt.append((me, ce))
                removed.append((max(cs, ms), min(ce, me)))
            cur = nxt
        for cs, ce in cur:
            if ce > cs:
                kept.append((cs, ce, name, strand))
    return kept, removed


def dustmasker_lowcplx(raw_bed, fasta, tmp):
    """Run dustmasker over merged target regions; return dict chrom->list[(s,e)]."""
    # merge regions
    regions = defaultdict(list)
    with open(raw_bed) as fh:
        for line in fh:
            f = line.split("\t")
            regions[f[0]].append((int(f[1]), int(f[2])))
    reg_bed = Path(tmp) / "regions.bed"
    with open(reg_bed, "w") as o:
        for chrom in sorted(regions):
            for s, e in merge_intervals(regions[chrom]):
                o.write(f"{chrom}\t{s}\t{e}\n")
    reg_fa = Path(tmp) / "regions.fa"
    run(["bedtools", "getfasta", "-fi", fasta, "-bed", str(reg_bed), "-fo", str(reg_fa)])
    dust = run(["dustmasker", "-in", str(reg_fa), "-outfmt", "acclist"])
    # acclist emits ">seqid\tstart\tend" per low-complexity range; seqid = chrom:start-end
    low = defaultdict(list)
    for line in dust.stdout.splitlines():
        f = line.lstrip(">").split("\t")
        if len(f) >= 3:
            m = re.match(r"(\S+):(\d+)-(\d+)", f[0])
            if m:
                chrom, off = m.group(1), int(m.group(2))
                low[chrom].append((off + int(f[1]), off + int(f[2]) + 1))
    return {c: merge_intervals(v) for c, v in low.items()}


def norm_chrom(c):
    """UCSC 'chr1'/'chrM' -> Ensembl '1'/'MT' to match the genome/GTF naming."""
    if c.startswith("chr"):
        c = c[3:]
    return "MT" if c in ("M", "MT") else c


def load_hard_mask(path, target_chroms):
    """Flexible loader: BED (0,1,2), UCSC rmsk (5,6,7), segdups/simpleRepeat (1,2,3).
    Normalizes UCSC chr-prefixed contigs to Ensembl naming."""
    if not path:
        return defaultdict(list)
    import gzip
    op = gzip.open if str(path).endswith(".gz") else open
    name = os.path.basename(path).lower()
    if "rmsk" in name:
        ci, si, ei = 5, 6, 7
    elif "superdup" in name or "segdup" in name or "simplerepeat" in name:
        ci, si, ei = 1, 2, 3
    else:
        ci, si, ei = 0, 1, 2
    out = defaultdict(list)
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            try:
                chrom = norm_chrom(f[ci])
                if chrom not in target_chroms:
                    continue
                out[chrom].append((int(f[si]), int(f[ei])))
            except (IndexError, ValueError):
                continue
    return {c: merge_intervals(v) for c, v in out.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--rmsk", default=None)
    ap.add_argument("--segdups", default=None)
    ap.add_argument("--simplerepeat", default=None, help="UCSC simpleRepeat track → STR repeat core")
    ap.add_argument("--relax-lowcplx-promoter", action="store_true", default=True)
    ap.add_argument("--str-skip-core", action="store_true", default=True)
    ap.add_argument("-o", "--output", required=True, help="targets.masked.bed")
    ap.add_argument("--lost", required=True, help="targets.lost_to_mask.bed")
    ap.add_argument("--summary", default=None)
    args = ap.parse_args()

    by_class = read_bed(args.raw)
    target_chroms = {c for cls in by_class.values() for c in cls}

    with tempfile.TemporaryDirectory() as tmp:
        low = dustmasker_lowcplx(args.raw, args.fasta, tmp)
    hard = load_hard_mask(args.rmsk, target_chroms)
    seg = load_hard_mask(args.segdups, target_chroms)
    simple = load_hard_mask(args.simplerepeat, target_chroms)

    def hard_for(chrom):
        return merge_intervals(hard.get(chrom, []) + seg.get(chrom, []))

    # policy: which masks to subtract per class
    #   hard = rmsk + segdups ; low = dustmasker low-complexity ; simple = tandem-repeat core
    POLICY = {
        "coding":     ("hard", "low"),
        "cnv":        ("hard", "low"),
        "promoter":   ("hard",),                       # relax low-complexity
        "regulatory": ("hard",),
        "str":        ("hard", "low", "simple"),       # strip repeat core -> unique flanks
        "anchor":     (),                              # minimal
        "extra":      (),
    }

    masked_out = open(args.output, "w")
    lost_out = open(args.lost, "w")
    str_core_out = open(re.sub(r"\.bed$", "", args.output) + ".str_core.bed", "w")

    bp_raw, bp_masked_kept, bp_removed = defaultdict(int), defaultdict(int), defaultdict(int)
    fully_lost = []
    per_target_kept = defaultdict(int)
    per_target_raw = defaultdict(int)

    for cls, chroms in by_class.items():
        which = POLICY.get(cls, ("hard", "low"))
        for chrom, ivs in chroms.items():
            for s, e, name, strand in ivs:
                bp_raw[cls] += e - s
                per_target_raw[name] += e - s
            masks = []
            if "hard" in which:
                masks += hard_for(chrom)
            if "low" in which:
                masks += low.get(chrom, [])
            if "simple" in which:
                masks += simple.get(chrom, [])
            kept, removed = subtract(ivs, masks)
            for s, e, name, strand in kept:
                masked_out.write(f"{chrom}\t{s}\t{e}\t{name}\t0\t{strand}\n")
                bp_masked_kept[cls] += e - s
                per_target_kept[name] += e - s
            for s, e in removed:
                lost_out.write(f"{chrom}\t{s}\t{e}\t{cls}_masked\t0\t.\n")
                bp_removed[cls] += e - s
                if cls == "str":
                    str_core_out.write(f"{chrom}\t{s}\t{e}\tstr_repeat_core\t0\t.\n")

    masked_out.close(); lost_out.close(); str_core_out.close()

    # flag targets fully lost to masking
    for name, raw_bp in per_target_raw.items():
        if per_target_kept.get(name, 0) == 0:
            fully_lost.append(name)

    summary = {
        "bp_raw_by_class": dict(bp_raw),
        "bp_kept_by_class": dict(bp_masked_kept),
        "bp_removed_by_class": dict(bp_removed),
        "frac_kept_by_class": {c: round(bp_masked_kept[c] / bp_raw[c], 4) for c in bp_raw if bp_raw[c]},
        "n_targets_fully_lost": len(fully_lost),
        "targets_fully_lost": fully_lost[:50],
        "hard_mask_used": bool(args.rmsk or args.segdups),
        "low_complexity_source": "dustmasker",
    }
    if args.summary:
        Path(args.summary).write_text(json.dumps(summary, indent=2))
    sys.stderr.write(f"[apply_masking] kept bp by class={dict(bp_masked_kept)} | "
                     f"fully_lost={len(fully_lost)}\n")
    if fully_lost:
        sys.stderr.write(f"[apply_masking] FLAG fully-masked targets (need boost/relax/accept): "
                         f"{fully_lost[:20]}\n")


if __name__ == "__main__":
    main()
