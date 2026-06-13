#!/usr/bin/env python3
"""
CaptureForge — Stage 8 : scaffold_synthesis.py

Wrap each unique bait sequence for the chosen synthesis format and validate.

Formats (--format):
  dna_pool     bait core only (minimal/no handle) — synthesized DNA oligo pool
  ivt_rna      5'-[PCR handle]-[T7 promoter]-[bait]-[PCR handle]-3' → in-vitro
               transcription to biotinylated RNA baits (amplifiable, low cost/capture)
  ready_to_use bait core only + targets (for a finished-bait vendor service)

Validates final length vs --max-len, rejects disallowed characters, and (ivt_rna)
checks T7 + both handles are present. Reverse-complements per --strand.

Outputs: synthesis.fasta (final oligos), synthesis.tsv, stage8_summary.json
"""
from __future__ import annotations
import argparse, csv, json, sys
from collections import defaultdict
from pathlib import Path

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

def revcomp(s):
    return s.translate(COMP)[::-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baits-tsv", required=True, help="baits.dedup.tsv (bait_id,sequence,copy_number,…)")
    ap.add_argument("--format", default="ivt_rna", choices=["dna_pool", "ivt_rna", "ready_to_use"])
    ap.add_argument("--t7", default="TAATACGACTCACTATAGGG")
    ap.add_argument("--handle5", default="ATCGCACCAGCGTGT")
    ap.add_argument("--handle3", default="CACTGCGGCTCCTCA")
    ap.add_argument("--max-len", type=int, default=170)
    ap.add_argument("--strand", default="plus", choices=["plus", "minus", "both"])
    ap.add_argument("--prefix", default="synthesis")
    ap.add_argument("--summary", default=None)
    args = ap.parse_args()

    fa = open(f"{args.prefix}.fasta", "w")
    tsv = open(f"{args.prefix}.tsv", "w")
    tsv.write("bait_id\tformat\tcore_len\tfinal_len\thas_t7\thas_handles\tover_max_len\tfinal_seq\n")

    lendist = defaultdict(int)
    n = over = bad = no_handle = 0
    for row in csv.DictReader(open(args.baits_tsv), delimiter="\t"):
        bait_id = row["bait_id"]; core = row["sequence"].upper()
        if args.strand == "minus":
            core = revcomp(core)
        if args.format == "ivt_rna":
            final = args.handle5 + args.t7 + core + args.handle3
        else:                                      # dna_pool / ready_to_use
            final = core
        has_t7 = args.t7 in final
        has_handles = (args.handle5 in final and args.handle3 in final) if args.format == "ivt_rna" else True
        over_max = len(final) > args.max_len
        bad_chars = any(c not in "ACGTU" for c in final)
        n += 1; lendist[len(final)] += 1
        if over_max: over += 1
        if bad_chars: bad += 1
        if args.format == "ivt_rna" and not (has_t7 and has_handles): no_handle += 1
        fa.write(f">{bait_id} fmt={args.format} len={len(final)} copies={row.get('copy_number','1')}\n{final}\n")
        tsv.write(f"{bait_id}\t{args.format}\t{len(core)}\t{len(final)}\t{has_t7}\t{has_handles}\t{over_max}\t{final}\n")
    fa.close(); tsv.close()

    summary = {
        "format": args.format, "n_baits": n,
        "final_length_distribution": {str(k): v for k, v in sorted(lendist.items())},
        "max_len": args.max_len, "n_over_max_len": over,
        "n_disallowed_chars": bad, "n_missing_handle": no_handle,
        "construct": (f"{args.handle5}+T7({args.t7})+bait+{args.handle3}"
                      if args.format == "ivt_rna" else "bait core only"),
        "pass": over == 0 and bad == 0 and no_handle == 0,
    }
    if args.summary:
        Path(args.summary).write_text(json.dumps(summary, indent=2))
    sys.stderr.write(f"[scaffold_synthesis] {n} baits fmt={args.format} | over_max={over} "
                     f"bad_chars={bad} missing_handle={no_handle} | lengths={dict(lendist)}\n")
    if over or bad or no_handle:
        sys.stderr.write("[scaffold_synthesis] ⚠ synthesis validation FLAGS present — review before ordering\n")


if __name__ == "__main__":
    main()
