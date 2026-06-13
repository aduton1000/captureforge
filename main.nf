#!/usr/bin/env nextflow
// ═══════════════════════════════════════════════════════════════════════════
//  CaptureForge — hybridization-capture probe (bait) design pipeline
//  main.nf  (entrypoint; full workflow wiring lands in task T10)
// ═══════════════════════════════════════════════════════════════════════════
nextflow.enable.dsl = 2

// NOTE: This is the T1 skeleton entrypoint. Subworkflows/modules for Stages 1–10
// are added in tasks T5–T10; this file currently validates inputs and prints the
// resolved plan so the scaffold is runnable and `-resume`-able from day one.

def helpMessage() {
    log.info """
    ╔══════════════════════════════════════════════════════════════════════╗
    ║  CaptureForge ${workflow.manifest.version}  —  bait design pipeline   ║
    ╚══════════════════════════════════════════════════════════════════════╝

    Usage:
      nextflow run main.nf -profile <test|mac_local|hpc_slurm>,<conda|docker|apptainer> \\
        --input <genes.txt|spec.csv|panel.xlsx> --genome_fasta <ref.fa> --gtf <ann.gtf>

    Key params (see params.example.yaml and nextflow.config for all):
      --input            gene list: plain-text symbols, CSV/TSV spec, or Excel workbook
      --genome_fasta     primary/no-alt reference FASTA   (required)
      --gtf              Ensembl GTF gene models          (required)
      --bait_len         bait length in nt                (default 120)
      --bait_format      dna_pool | ivt_rna | ready_to_use
      --outdir           results directory                (default results/)
    """.stripIndent()
}

workflow {
    if (params.help) { helpMessage(); return }

    // Minimal input validation so the scaffold fails loudly & early.
    if (!params.input)        { exit 1, "ERROR: --input gene list is required" }
    if (!params.genome_fasta) { exit 1, "ERROR: --genome_fasta is required" }

    log.info """
    ─────────────────────────────────────────────────────────────────────────
    CaptureForge ${workflow.manifest.version}   (scaffold entrypoint)
    ─────────────────────────────────────────────────────────────────────────
    panel        : ${params.panel_name}
    input        : ${params.input}
    genome_build : ${params.genome_build}
    genome_fasta : ${params.genome_fasta}
    gtf          : ${params.gtf}
    bait_len     : ${params.bait_len} nt
    bait_format  : ${params.bait_format}
    outdir       : ${params.outdir}
    profile      : ${workflow.profile}
    ─────────────────────────────────────────────────────────────────────────
    Stage modules (T5–T10) not yet wired. This entrypoint will invoke
    subworkflows/captureforge.nf once implemented.
    """.stripIndent()
}
