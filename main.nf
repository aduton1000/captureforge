#!/usr/bin/env nextflow
// ═══════════════════════════════════════════════════════════════════════════
//  CaptureForge — hybridization-capture probe (bait) design pipeline
//  main.nf  (entrypoint → subworkflows/captureforge.nf)
// ═══════════════════════════════════════════════════════════════════════════
nextflow.enable.dsl = 2

include { CAPTUREFORGE } from './subworkflows/captureforge.nf'

def helpMessage() {
    log.info """
    ╔══════════════════════════════════════════════════════════════════════╗
    ║  CaptureForge ${workflow.manifest.version}  —  bait design pipeline   ║
    ╚══════════════════════════════════════════════════════════════════════╝
    Usage:
      nextflow run main.nf -profile <test|local|hpc_slurm>,<conda|docker|apptainer> \\
        --input <genes.txt|spec.csv|panel.xlsx> --genome_fasta <noalt.fa> --gtf <ann.gtf> \\
        --specificity_db <noalt_index_prefix>

    Required: --input  --genome_fasta  --gtf  --specificity_db
    See params.example.yaml / nextflow.config for all parameters.
    """.stripIndent()
}

workflow {
    if (params.help) { helpMessage(); return }
    if (!params.input)         { exit 1, "ERROR: --input gene list is required" }
    if (!params.genome_fasta)  { exit 1, "ERROR: --genome_fasta is required" }
    if (!params.gtf)           { exit 1, "ERROR: --gtf is required" }
    if (!params.specificity_db) {
        log.warn "No --specificity_db set: Stage 4/6 will fail. Build one from the NO-ALT genome:\n" +
                 "  bin/build_specificity_index.sh --fasta ${params.genome_fasta} --aligner ${params.specificity_aligner} --out refs/idx"
    }

    log.info """
    ─────────────────────────────────────────────────────────────────────────
    CaptureForge ${workflow.manifest.version}   profile=${workflow.profile}
    panel=${params.panel_name}  build=${params.genome_build}  bait_len=${params.bait_len}  fmt=${params.bait_format}
    input=${params.input}
    genome=${params.genome_fasta}
    gtf=${params.gtf}   specificity_db=${params.specificity_db}
    outdir=${params.outdir}
    ─────────────────────────────────────────────────────────────────────────
    """.stripIndent()

    gene_list = file(params.input,        checkIfExists: true)
    gtf       = file(params.gtf,          checkIfExists: true)
    genome    = file(params.genome_fasta, checkIfExists: true)

    CAPTUREFORGE(gene_list, gtf, genome)
}
