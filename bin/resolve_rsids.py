#!/usr/bin/env python3
"""
CaptureForge — resolve_rsids.py

Resolve rsIDs to GRCh38 coordinates via Ensembl Variation (REST), emitting an
rsID BED that define_targets.py consumes through --rsid-bed. This is the
lightweight path (no 20 GB dbSNP VCF) and uses an Ensembl-listed source; for a
pinned production run on the HPC, fetch the dbSNP VCF instead
(bin/fetch_references.sh dbsnp) and pass --dbsnp to define_targets.

Usage:
  resolve_rsids.py --from-targets targets.tsv -o anchors.rsid.bed \
     [--manifest refs/refs.provenance.jsonl] [--rest https://rest.ensembl.org]
  resolve_rsids.py --rsids rs123,rs456 -o anchors.rsid.bed
"""
from __future__ import annotations
import argparse, csv, hashlib, json, re, sys, time, urllib.request
from pathlib import Path

RSID_RE = re.compile(r"\brs\d+\b", re.I)


def collect_rsids(args):
    ids = []
    if args.rsids:
        ids += RSID_RE.findall(args.rsids)
    if args.from_targets:
        for row in csv.DictReader(open(args.from_targets), delimiter="\t"):
            ids += RSID_RE.findall(row.get("anchors", "") or "")
    # dedupe, preserve order
    seen, out = set(), []
    for i in ids:
        k = i.lower()
        if k not in seen:
            seen.add(k); out.append(i)
    return out


def rest_post(rest, ids, species="homo_sapiens", chunk=200):
    """POST /variation/{species}; returns dict rsid->record."""
    out = {}
    for i in range(0, len(ids), chunk):
        batch = ids[i:i + chunk]
        body = json.dumps({"ids": batch}).encode()
        req = urllib.request.Request(
            f"{rest}/variation/{species}",
            data=body, method="POST",
            headers={"Content-Type": "application/json", "Accept": "application/json"})
        for attempt in range(4):
            try:
                with urllib.request.urlopen(req, timeout=60) as r:
                    out.update(json.load(r))
                break
            except Exception as e:
                if attempt == 3:
                    sys.stderr.write(f"[resolve_rsids] REST error: {e}\n")
                else:
                    time.sleep(2 * (attempt + 1))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from-targets")
    ap.add_argument("--rsids")
    ap.add_argument("--rest", default="https://rest.ensembl.org")
    ap.add_argument("--assembly", default="GRCh38")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--manifest", help="append a provenance record (JSONL) here")
    args = ap.parse_args()

    ids = collect_rsids(args)
    if not ids:
        sys.exit("ERROR: no rsIDs found")
    sys.stderr.write(f"[resolve_rsids] resolving {len(ids)} rsIDs via {args.rest} …\n")
    data = rest_post(args.rest, ids)

    rows, unresolved = [], []
    for rs in ids:
        rec = data.get(rs)
        mapping = None
        if rec:
            for m in rec.get("mappings", []):
                if m.get("assembly_name") == args.assembly:
                    mapping = m; break
        if not mapping:
            unresolved.append(rs); continue
        chrom = str(mapping["seq_region_name"])
        start = int(mapping["start"])            # 1-based
        pos0 = start - 1
        rows.append((chrom, pos0, start, rs.lower()))

    rows.sort(key=lambda r: (r[0], r[1]))
    with open(args.output, "w") as o:
        for chrom, s, e, rs in rows:
            o.write(f"{chrom}\t{s}\t{e}\t{rs}\n")

    sys.stderr.write(f"[resolve_rsids] resolved {len(rows)}/{len(ids)} → {args.output}"
                     + (f" | unresolved: {unresolved}" if unresolved else "") + "\n")

    if args.manifest:
        md5 = hashlib.md5(Path(args.output).read_bytes()).hexdigest()
        rec = {"name": "dbsnp_anchors_ensembl_rest", "build": args.assembly,
               "release": "Ensembl Variation (REST, current)",
               "source": f"{args.rest}/variation/homo_sapiens",
               "path": args.output, "n_resolved": len(rows), "n_requested": len(ids),
               "md5": md5}
        with open(args.manifest, "a") as m:
            m.write(json.dumps(rec) + "\n")
        sys.stderr.write(f"[resolve_rsids] provenance appended → {args.manifest}\n")


if __name__ == "__main__":
    main()
