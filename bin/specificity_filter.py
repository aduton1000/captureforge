#!/usr/bin/env python3
"""
CaptureForge — Stage 4 : specificity_filter.py

Detect off-target baits by aligning every candidate to the whole genome and
counting secondary hits above identity/length thresholds. A bait that maps
uniquely (only its own locus) passes; one with qualifying off-target hits is
rejected (or flagged, depending on tolerance).

Backends (--aligner):
  blast    blastn megablast vs a BLAST nucleotide DB  (no index build; used locally)
  bowtie2  bowtie2 -a vs a prebuilt index             (HPC default)
  bwa      bwa mem -a vs a prebuilt index

Self-hit handling is naming-agnostic: each bait has exactly one perfect,
full-length hit (its source locus); off_target_count = (qualifying hits) - 1.

Outputs: candidates.spec.tsv (per-bait verdict), candidates.pass.bed, summary JSON.
"""
from __future__ import annotations
import argparse, json, os, re, subprocess, sys
from collections import defaultdict
from pathlib import Path


# Contigs that must NOT be in a specificity reference: alt haplotypes, decoys,
# HLA/MHC alternates. Unplaced/unlocalized scaffolds (KI*/GL*) ARE part of the
# no-alt analysis set and are allowed.
ALT_RE = re.compile(r"(_alt\b|HSCHR|_MHC_|\bMHC_|HLA[-_]|_decoy|hs38d1|chrUn_.*decoy|EBV)", re.I)


def db_seqnames(db, aligner):
    """Best-effort enumeration of the specificity reference's sequence names."""
    try:
        if aligner == "blast":
            out = subprocess.run(["blastdbcmd", "-db", db, "-entry", "all", "-outfmt", "%t"],
                                 capture_output=True, text=True, timeout=120)
            return [ln.split()[0] for ln in out.stdout.splitlines() if ln.strip()]
        if aligner == "bowtie2":
            out = subprocess.run(["bowtie2-inspect", "-n", db], capture_output=True, text=True, timeout=120)
            return [ln.strip() for ln in out.stdout.splitlines() if ln.strip()]
        if aligner == "bwa":                                  # names live in <db>.ann
            ann = db + ".ann"
            if os.path.exists(ann):
                names = []
                with open(ann) as fh:
                    next(fh, None)
                    for i, ln in enumerate(fh):
                        if i % 2 == 0 and ln.split():
                            names.append(ln.split()[0])
                return names
    except Exception as e:
        sys.stderr.write(f"[specificity] WARN could not enumerate db contigs ({e}); ref-invariant check skipped\n")
    return []


def assert_noalt_reference(db, aligner, genome_fai=None):
    """FAIL LOUD if the specificity reference contains alt/decoy/HLA contigs, or
    if it disagrees with the design genome. This prevents the classic trap of
    pointing CaptureForge at a convenient toplevel/with-alt BLAST db, which makes
    MHC-region genes (e.g. GENE) look non-specific vs their own alt haplotypes."""
    import os
    names = db_seqnames(db, aligner)
    if not names:
        return
    alts = sorted({n for n in names if ALT_RE.search(n)})
    if alts:
        sys.stderr.write(
            "\nERROR: specificity reference contains ALT/decoy/HLA contigs — specificity MUST run against "
            "the NO-ALT primary assembly (the same reference as --genome_fasta). Offending contigs:\n"
            + "".join(f"   - {a}\n" for a in alts[:20])
            + (f"   … (+{len(alts)-20} more)\n" if len(alts) > 20 else "")
            + "Rebuild from the no-alt genome:\n"
              f"   bin/build_specificity_index.sh --fasta <primary_no_alt.fa> --aligner {aligner} --out <prefix>\n")
        sys.exit(3)
    if genome_fai and os.path.exists(genome_fai):
        gnames = {l.split('\t')[0] for l in open(genome_fai)}
        dnames = set(names)
        extra = dnames - gnames
        if extra:
            sys.stderr.write(
                f"ERROR: specificity db and --genome_fasta are DIFFERENT references "
                f"({len(extra)} db contigs absent from the genome, e.g. {sorted(extra)[:5]}). "
                f"They must be the same no-alt reference.\n")
            sys.exit(3)
    sys.stderr.write(f"[specificity] reference invariant OK: {len(names)} contigs, no alt/decoy/HLA.\n")


def read_bed(path):
    rows = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 4:
            rows[f[3]] = line
    return rows


def _clean(name):
    return re.sub(r"\([+-]\)$", "", name)   # defensive: drop any strand suffix


def qlen_from_fasta(path):
    lens, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = _clean(line[1:].strip().split()[0])
            lens[name] = 0
        elif name:
            lens[name] += len(line.strip())
    return lens


def run_blast(fasta, db, threads, min_id):
    cmd = ["blastn", "-task", "megablast", "-query", fasta, "-db", db,
           "-outfmt", "6 qseqid sseqid pident length qlen sstart send",
           "-perc_identity", str(int(min_id * 100 - 10)),  # pre-filter a bit below threshold
           "-max_target_seqs", "100", "-dust", "no", "-num_threads", str(threads)]
    sys.stderr.write(f"[specificity] blastn megablast (threads={threads}) …\n")
    p = subprocess.run(cmd, check=True, text=True, capture_output=True)
    hits = defaultdict(list)  # qseqid -> list[(pident, length, qlen)]
    for line in p.stdout.splitlines():
        f = line.split("\t")
        hits[_clean(f[0])].append((float(f[2]), int(f[3]), int(f[4])))
    return hits


def run_sam(fasta, index, threads, aligner):
    if aligner == "bowtie2":
        cmd = ["bowtie2", "-x", index, "-f", "-U", fasta, "-a", "--quiet",
               "-p", str(threads), "--no-hd"]
    else:  # bwa mem
        cmd = ["bash", "-c", f"bwa mem -a -t {threads} {index} {fasta}"]
    sys.stderr.write(f"[specificity] {aligner} all-hits alignment …\n")
    p = subprocess.run(cmd, check=True, text=True, capture_output=True)
    hits = defaultdict(list)
    for line in p.stdout.splitlines():
        if line.startswith("@") or not line.strip():
            continue
        f = line.split("\t")
        qid, flag, cigar, seq = _clean(f[0]), int(f[1]), f[5], f[9]
        if flag & 4:
            continue
        qlen = len(seq) if seq != "*" else sum(int(n) for n, op in re.findall(r"(\d+)([MIS=X])", cigar))
        m = sum(int(n) for n, op in re.findall(r"(\d+)([M=X])", cigar))
        nm = next((int(t.split(":")[-1]) for t in f[11:] if t.startswith("NM:i:")), 0)
        ident = 100.0 * (m - nm) / m if m else 0.0
        hits[qid].append((ident, m, qlen))
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True, help="candidates.fasta")
    ap.add_argument("--candidates-bed", required=True)
    ap.add_argument("--aligner", default="blast", choices=["blast", "bowtie2", "bwa"])
    ap.add_argument("--db", help="BLAST db prefix (blast) or aligner index prefix")
    ap.add_argument("--min-identity", type=float, default=0.90)
    ap.add_argument("--min-lenfrac", type=float, default=0.50)
    ap.add_argument("--max-offtarget", type=int, default=1,
                    help="tolerated off-target hits; reject if off_target_count > this")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--genome-fai", default=None, help="design genome .fai for same-reference cross-check")
    ap.add_argument("--no-ref-check", action="store_true", help="skip the no-alt reference invariant (not recommended)")
    ap.add_argument("-o", "--output", required=True, help="candidates.pass.bed")
    ap.add_argument("--verdicts", required=True, help="candidates.spec.tsv")
    ap.add_argument("--summary", default=None)
    args = ap.parse_args()

    if not args.no_ref_check:
        assert_noalt_reference(args.db, args.aligner, args.genome_fai)

    bed = read_bed(args.candidates_bed)
    qlens = qlen_from_fasta(args.candidates)

    if args.aligner == "blast":
        hits = run_blast(args.candidates, args.db, args.threads, args.min_identity)
    else:
        hits = run_sam(args.candidates, args.db, args.threads, args.aligner)

    min_id_pct = args.min_identity * 100
    out = open(args.output, "w")
    vout = open(args.verdicts, "w")
    vout.write("bait_id\tclass\toff_target_count\tmax_offhit_id\tverdict\n")

    by_class_total = defaultdict(int)
    by_class_reject = defaultdict(int)
    offdist = defaultdict(int)
    n_pass = 0

    for bait_id, bedline in bed.items():
        cls = bait_id.split("|")[1] if "|" in bait_id else "coding"
        by_class_total[cls] += 1
        qlen = qlens.get(bait_id, 120)
        qhits = hits.get(bait_id, [])
        qualifying = [h for h in qhits if h[0] >= min_id_pct and h[1] >= args.min_lenfrac * qlen]
        # exclude one perfect/best self-hit
        off_count = max(0, len(qualifying) - 1)
        off_ids = sorted([h[0] for h in qualifying], reverse=True)[1:]   # identities of off-targets
        max_off = round(off_ids[0], 1) if off_ids else 0.0
        offdist[off_count] += 1
        verdict = "PASS" if off_count <= args.max_offtarget else "REJECT"
        vout.write(f"{bait_id}\t{cls}\t{off_count}\t{max_off}\t{verdict}\n")
        if verdict == "PASS":
            out.write(bedline)
            n_pass += 1
        else:
            by_class_reject[cls] += 1

    out.close(); vout.close()

    summary = {
        "aligner": args.aligner,
        "thresholds": {"min_identity": args.min_identity, "min_lenfrac": args.min_lenfrac,
                       "max_offtarget": args.max_offtarget},
        "n_candidates": sum(by_class_total.values()),
        "n_pass": n_pass,
        "n_reject": sum(by_class_reject.values()),
        "reject_by_class": dict(by_class_reject),
        "total_by_class": dict(by_class_total),
        "offtarget_count_distribution": {str(k): v for k, v in sorted(offdist.items())},
    }
    if args.summary:
        Path(args.summary).write_text(json.dumps(summary, indent=2))
    sys.stderr.write(f"[specificity] pass={n_pass}/{sum(by_class_total.values())} "
                     f"reject_by_class={dict(by_class_reject)}\n")


if __name__ == "__main__":
    main()
