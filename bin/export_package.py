#!/usr/bin/env python3
"""
CaptureForge — Stage 10 : export_package.py

Assemble the order-ready package: baits.fasta (final synthesis constructs),
baits.csv (full per-bait metadata), a generic order table, a MANIFEST of md5
checksums, and a run-provenance file (tool versions, reference releases, params,
git commit). CaptureForge does NOT contact any vendor — map baits.csv to a
specific synthesizer template before ordering (see assets/vendor_templates).
"""
from __future__ import annotations
import argparse, csv, hashlib, json, os, re, subprocess, sys
from pathlib import Path


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def tm_calc(s):
    try:
        import primer3
        return round(float(primer3.calc_tm(s)), 1)
    except Exception:
        return ""


def tool_version(cmd):
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
        return (out.stdout + out.stderr).strip().splitlines()[0]
    except Exception:
        return "n/a"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baits-tsv", required=True, help="baits.dedup.tsv")
    ap.add_argument("--baits-bed", required=True, help="baits.dedup.bed (strand/coords)")
    ap.add_argument("--synthesis-tsv", required=True, help="synthesis.tsv (final constructs)")
    ap.add_argument("--spec-tsv", default=None, help="candidates.spec.tsv (off-target flag join)")
    ap.add_argument("--format", default="ivt_rna")
    ap.add_argument("--genome", default=None, help="design genome FASTA path (recorded in provenance)")
    ap.add_argument("--specificity-db", default=None, help="specificity reference/db (recorded in provenance)")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--manifest-files", default="", help="comma-sep extra files to checksum")
    ap.add_argument("--params-json", default=None, help="JSON of run params to record")
    ap.add_argument("--refs-provenance", default=None, help="refs.provenance.json to fold in")
    ap.add_argument("--summary", default=None)
    args = ap.parse_args()

    od = Path(args.outdir); od.mkdir(parents=True, exist_ok=True)

    # strand + coords per bait_id
    strand, coord = {}, {}
    for l in open(args.baits_bed):
        f = l.rstrip("\n").split("\t")
        strand[f[3]] = f[5] if len(f) > 5 else "+"
        coord[f[3]] = f"{f[0]}:{f[1]}-{f[2]}"

    # final construct per bait_id
    final_seq = {}
    for row in csv.DictReader(open(args.synthesis_tsv), delimiter="\t"):
        final_seq[row["bait_id"]] = row["final_seq"]

    # off-target flag: any source candidate had a (tolerated) off-target hit
    off_count = {}
    if args.spec_tsv and os.path.exists(args.spec_tsv):
        for row in csv.DictReader(open(args.spec_tsv), delimiter="\t"):
            off_count[row["bait_id"]] = int(row["off_target_count"])

    baits_csv = open(od / "baits.csv", "w", newline="")
    w = csv.writer(baits_csv)
    w.writerow(["bait_id", "sequence", "copy_number", "gene", "class", "source_region",
                "strand", "GC", "Tm", "off_target_flag"])
    baits_fa = open(od / "baits.fasta", "w")
    order = open(od / "order_generic.tsv", "w")
    order.write("Pool_Name\tProbe_ID\tSequence\tReplicates\n")

    n = 0
    for row in csv.DictReader(open(args.baits_tsv), delimiter="\t"):
        bid = row["bait_id"]; core = row["sequence"]
        fseq = final_seq.get(bid, core)
        # off-target flag via source candidate ids
        flag = False
        for src in row.get("source_regions", "").split(";"):
            name = src.split(":")[-1] if ":" in src else src
            if off_count.get(name, 0) >= 1:
                flag = True; break
        w.writerow([bid, fseq, row.get("copy_number", "1"), row.get("genes", ""),
                    row.get("classes", ""), (row.get("source_regions", "").split(";")[0]),
                    strand.get(bid, "+"), row.get("gc", ""), tm_calc(core), int(flag)])
        baits_fa.write(f">{bid} copies={row.get('copy_number','1')} fmt={args.format}\n{fseq}\n")
        order.write(f"CaptureForge_pool\t{bid}\t{fseq}\t{row.get('copy_number','1')}\n")
        n += 1
    baits_csv.close(); baits_fa.close(); order.close()

    # ── MANIFEST (md5 of key outputs) ──
    files = ["baits.fasta", "baits.csv", "order_generic.tsv"]
    files += [f for f in args.manifest_files.split(",") if f]
    manifest = od / "MANIFEST"
    with open(manifest, "w") as m:
        for f in files:
            p = Path(f)
            if not p.is_absolute() and (od / f).exists():
                p = od / f
            if p.exists():
                m.write(f"{md5(p)}  {p.name}  ({p.stat().st_size} bytes)\n")

    # ── provenance ──
    prov = {
        "tool_versions": {
            "bedtools": tool_version(["bedtools", "--version"]),
            "samtools": tool_version(["samtools", "--version"]),
            "blastn": tool_version(["blastn", "-version"]),
            "bowtie2": tool_version(["bowtie2", "--version"]),
            "dustmasker": tool_version(["dustmasker", "-version"]),
            "python": sys.version.split()[0],
        },
        "git_commit": tool_version(["git", "rev-parse", "HEAD"]),
        "bait_format": args.format,
        "n_unique_baits": n,
        "design_genome": args.genome,
        "specificity_reference": args.specificity_db,
        "reference_invariant": "specificity ran against the no-alt primary assembly == design genome "
                               "(enforced by specificity_filter.py assert_noalt_reference)",
    }
    if args.params_json and os.path.exists(args.params_json):
        try:
            prov["params"] = json.loads(Path(args.params_json).read_text())
        except Exception:
            pass
    if args.refs_provenance and os.path.exists(args.refs_provenance):
        try:
            prov["references"] = json.loads(Path(args.refs_provenance).read_text())
        except Exception:
            pass
    (od / "provenance.json").write_text(json.dumps(prov, indent=2))

    summary = {"n_baits_exported": n, "format": args.format,
               "files": [f for f in files] + ["MANIFEST", "provenance.json"]}
    if args.summary:
        Path(args.summary).write_text(json.dumps(summary, indent=2))
    sys.stderr.write(f"[export_package] exported {n} baits → baits.fasta / baits.csv / order_generic.tsv "
                     f"+ MANIFEST + provenance.json\n")


if __name__ == "__main__":
    main()
