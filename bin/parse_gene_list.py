#!/usr/bin/env python3
"""
CaptureForge — Stage 1a : parse_gene_list.py

Normalize ANY accepted gene-list on-ramp (plain text / CSV-TSV spec / Excel
workbook) into the internal target table consumed by the rest of the pipeline.

The logic that turns genes into labelled, design-classed targets is part of the
pipeline itself: a brand-new user with nothing but a list of symbols can run it.

Usage:
    parse_gene_list.py --input <file> --gtf <Ensembl.gtf[.gz]> \
        [--excel-sheet NAME] [--default-classes coding,promoter] \
        [--anchor-flank 150] [--allow-unresolved] \
        -o targets.tsv --summary parse_summary.json

Outputs the internal table (see assets/gene_list_schema.md):
    gene  ensembl_id  chrom  strand  gene_start  gene_end
    design_classes  anchors  extra_regions  burden_group
"""
from __future__ import annotations
import argparse, csv, gzip, json, re, sys
from pathlib import Path

RSID_RE = re.compile(r"\brs\d+\b", re.IGNORECASE)
TRUE_TOKENS = {"1", "yes", "y", "true", "t", "coding"}

# ── Variant-class text → design classes ────────────────────────────────────
CLASS_RULES = [
    ("regulatory", re.compile(r"enhancer|distal|regulatory", re.I)),
    ("promoter",   re.compile(r"promoter|\breg\b|regulat|enhancer|haplotype|upstream", re.I)),
    ("cnv",        re.compile(r"\bcnv\b|copy.?number|structural|deletion|duplication", re.I)),
    ("str",        re.compile(r"\bstr\b|repeat|tandem|microsat|expansion", re.I)),
]

def map_design_classes(text: str, base: set[str] | None = None) -> list[str]:
    classes = set(base) if base else set()
    classes.add("coding")                       # coding is always implied
    if text:
        for cls, rx in CLASS_RULES:
            if rx.search(text):
                classes.add(cls)
    # deterministic order
    order = ["coding", "promoter", "regulatory", "cnv", "str", "anchor"]
    return [c for c in order if c in classes]


# ── GTF gene index (symbol/alias → ensembl id + coordinates) ───────────────
def _openmaybe(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)

ATTR_NAME = re.compile(r'gene_name "([^"]+)"')
ATTR_ID   = re.compile(r'gene_id "([^"]+)"')

def load_gtf_genes(gtf: str):
    """Return (name2rec, lc2name) where rec = dict(ensembl_id,chrom,strand,start,end)."""
    name2rec: dict[str, dict] = {}
    lc2name: dict[str, str] = {}
    with _openmaybe(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            mn = ATTR_NAME.search(f[8]); mi = ATTR_ID.search(f[8])
            if not mn or not mi:
                continue
            name = mn.group(1)
            rec = {"ensembl_id": mi.group(1), "chrom": f[0], "strand": f[6],
                   "start": int(f[3]), "end": int(f[4])}
            name2rec[name] = rec
            lc2name.setdefault(name.lower(), name)
    return name2rec, lc2name


def resolve_symbol(raw: str, name2rec: dict, lc2name: dict):
    """Resolve 'A (B)' / aliased symbol to a valid GTF gene_name. Returns (name, rec) or (None, None)."""
    raw = (raw or "").strip()
    if not raw:
        return None, None
    cands: list[str] = []
    m = re.match(r"^(.*?)\s*\(([^)]+)\)\s*$", raw)
    if m:
        cands += [m.group(1).strip(), m.group(2).strip()]
    else:
        cands.append(raw)
    # also split on slash / comma (e.g. "A/B")
    cands = [c.strip() for part in cands for c in re.split(r"[/,]", part) if c.strip()]
    for c in cands:                              # exact match first
        if c in name2rec:
            return c, name2rec[c]
    for c in cands:                              # case-insensitive fallback
        hit = lc2name.get(c.lower())
        if hit:
            return hit, name2rec[hit]
    return None, None


# ── Input readers (auto-detected) ──────────────────────────────────────────
def detect_format(path: str) -> str:
    ext = Path(path).suffix.lower()
    if ext in (".xlsx", ".xlsm"):
        return "excel"
    if ext in (".csv", ".tsv"):
        return "delim"
    return "plain"

def read_plain(path: str):
    rows = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s and not s.startswith("#"):
                rows.append({"gene": s})
    return rows

def read_delim(path: str):
    delim = "\t" if path.lower().endswith(".tsv") else ","
    with open(path, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter=delim)
        rdr.fieldnames = [(h or "").strip().lower() for h in (rdr.fieldnames or [])]
        if "gene" not in rdr.fieldnames:
            sys.exit(f"ERROR: CSV/TSV spec must have a 'gene' column; found {rdr.fieldnames}")
        return [{(k or "").strip().lower(): (v.strip() if isinstance(v, str) else v)
                 for k, v in row.items()} for row in rdr]

def read_excel(path: str, sheet: str):
    try:
        from openpyxl import load_workbook
    except ImportError:
        sys.exit("ERROR: reading .xlsx requires openpyxl (pip install openpyxl)")
    wb = load_workbook(path, read_only=True, data_only=True)
    ws = wb[sheet] if sheet in wb.sheetnames else wb[wb.sheetnames[0]]
    grid = [[("" if c is None else str(c)).strip() for c in row] for row in ws.iter_rows(values_only=True)]
    # find header row = the one containing a cell exactly 'Gene'
    hdr_idx = next((i for i, r in enumerate(grid) if any(c.lower() == "gene" for c in r)), None)
    if hdr_idx is None:
        sys.exit(f"ERROR: no 'Gene' header found in sheet '{ws.title}'")
    header = [c.strip() for c in grid[hdr_idx]]
    col = {h.lower(): j for j, h in enumerate(header) if h}
    rows = []
    for r in grid[hdr_idx + 1:]:
        def get(name):
            j = col.get(name.lower())
            return r[j] if (j is not None and j < len(r)) else ""
        gene = get("Gene")
        if not gene:
            continue
        rows.append({
            "gene": gene,
            "variant class": get("Variant class"),
            "key variants to capture": get("Key variants to capture"),
            "locus": get("Locus"),
            "burden group": get("Burden group"),
        })
    return rows


def truthy(v) -> bool:
    return str(v).strip().lower() in TRUE_TOKENS if v not in (None, "") else False


def main():
    ap = argparse.ArgumentParser(description="Normalize a gene list into the CaptureForge target table.")
    ap.add_argument("--input", required=True)
    ap.add_argument("--gtf", required=True, help="Ensembl GTF (.gtf or .gtf.gz) for symbol→Ensembl resolution")
    ap.add_argument("--excel-sheet", default="Candidate Panel (expanded)")
    ap.add_argument("--default-classes", default="coding,promoter")
    ap.add_argument("--anchor-flank", type=int, default=150)
    ap.add_argument("--allow-unresolved", action="store_true")
    ap.add_argument("-o", "--output", default="-")
    ap.add_argument("--summary", default=None)
    args = ap.parse_args()

    fmt = detect_format(args.input)
    raw_rows = {"plain": read_plain, "delim": read_delim,
                "excel": lambda p: read_excel(p, args.excel_sheet)}[fmt](args.input)
    if not raw_rows:
        sys.exit(f"ERROR: no genes found in {args.input}")

    sys.stderr.write(f"[parse_gene_list] format={fmt}  rows={len(raw_rows)}  loading GTF…\n")
    name2rec, lc2name = load_gtf_genes(args.gtf)
    sys.stderr.write(f"[parse_gene_list] GTF genes indexed: {len(name2rec)}\n")

    default_base = {c.strip() for c in args.default_classes.split(",") if c.strip()}
    out_rows, unresolved, seen = [], [], set()

    for row in raw_rows:
        raw_gene = row.get("gene", "")
        name, rec = resolve_symbol(raw_gene, name2rec, lc2name)
        if not rec:
            unresolved.append(raw_gene)
            continue
        if name in seen:                          # de-duplicate
            continue
        seen.add(name)

        # design classes
        #   plain text : no per-gene info → apply the default_classes (e.g. coding+promoter)
        #   excel      : derive from the curated 'Variant class' text (coding always implied)
        #   csv/tsv    : coding always + whichever optional class columns are truthy
        if fmt == "excel":
            classes = set(map_design_classes(row.get("variant class", ""), base=None))
        elif fmt == "delim":
            classes = set(map_design_classes("", base=None))
            for c in ("promoter", "regulatory", "cnv", "str"):
                if truthy(row.get(c)):
                    classes.add(c)
        else:  # plain
            classes = set(default_base) | {"coding"}

        # anchors (rsIDs)
        anchor_text = " ".join(str(row.get(k, "")) for k in
                               ("key variants to capture", "anchors", "anchor"))
        anchors = sorted(set(m.group(0).lower() for m in RSID_RE.finditer(anchor_text)))
        if anchors:
            classes.add("anchor")

        order = ["coding", "promoter", "regulatory", "cnv", "str", "anchor"]
        classes = [c for c in order if c in classes]

        out_rows.append({
            "gene": name,
            "ensembl_id": rec["ensembl_id"],
            "chrom": rec["chrom"],
            "strand": rec["strand"],
            "gene_start": rec["start"],
            "gene_end": rec["end"],
            "design_classes": "|".join(classes),
            "anchors": ",".join(anchors),
            "extra_regions": str(row.get("extra_regions", "") or ""),
            "burden_group": str(row.get("burden group", "") or ""),
        })

    if unresolved and not args.allow_unresolved:
        sys.stderr.write("\nERROR: unresolvable gene symbols (not found in GTF gene_name):\n")
        for g in unresolved:
            sys.stderr.write(f"   - {g}\n")
        sys.stderr.write("Fix the symbols, or re-run with --allow-unresolved to skip them.\n")
        sys.exit(2)
    if unresolved:
        sys.stderr.write(f"[parse_gene_list] WARNING skipped {len(unresolved)} unresolved: {unresolved}\n")

    cols = ["gene", "ensembl_id", "chrom", "strand", "gene_start", "gene_end",
            "design_classes", "anchors", "extra_regions", "burden_group"]
    out = sys.stdout if args.output == "-" else open(args.output, "w", newline="")
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for r in out_rows:
        w.writerow(r)
    if out is not sys.stdout:
        out.close()

    # class tally + summary
    tally: dict[str, int] = {}
    for r in out_rows:
        for c in r["design_classes"].split("|"):
            tally[c] = tally.get(c, 0) + 1
    summary = {
        "input": args.input, "format": fmt,
        "n_input": len(raw_rows), "n_resolved": len(out_rows),
        "n_unresolved": len(unresolved), "unresolved": unresolved,
        "class_counts": tally,
        "n_with_anchors": sum(1 for r in out_rows if r["anchors"]),
    }
    if args.summary:
        Path(args.summary).write_text(json.dumps(summary, indent=2))
    sys.stderr.write(
        f"[parse_gene_list] resolved {len(out_rows)}/{len(raw_rows)} genes | "
        f"classes={tally} | anchors_on={summary['n_with_anchors']}\n")


if __name__ == "__main__":
    main()
