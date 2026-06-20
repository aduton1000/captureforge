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
      nextflow run main.nf -profile <test|local|hpc_slurm>,<conda|docker|apptainer|singularity> \\
        --input <genes.txt|spec.csv|panel.xlsx> --genome_fasta <noalt.fa> --gtf <ann.gtf> \\
        --specificity_db <noalt_index_prefix>

    Required: --input  --genome_fasta  --gtf  --specificity_db
    See params.example.yaml / nextflow.config for all parameters.
    """.stripIndent()
}

workflow {
    if (params.help) { helpMessage(); return }
    // Guard: an enabled container engine with no image silently runs tasks on the
    // HOST (host tools may be incomplete — e.g. no matplotlib). Fail loudly instead.
    if (workflow.containerEngine && !params.container_image) {
        exit 1, "ERROR: container engine '${workflow.containerEngine}' is enabled but --container_image is unset.\n" +
                "  Set it to the image so tasks run inside the container, e.g.\n" +
                "    singularity/apptainer:  --container_image /path/to/captureforge.sif\n" +
                "    docker:                 --container_image captureforge:0.1.0\n" +
                "  (or omit the engine profile to run natively on the host)."
    }
    // Guard: under the SLURM executor (profile hpc_slurm), submitting with no
    // partition relies on a default that may not exist on this cluster (e.g. the
    // old 'normal' default is rejected where the queue is 'global'). Fail loudly
    // rather than letting sbatch reject every job. Only fires for hpc_slurm — local
    // and test have no partition. (Short-circuits before reading slurm_partition.)
    if (workflow.profile.tokenize(',').contains('hpc_slurm') && !params.slurm_partition) {
        exit 1, "ERROR: SLURM executor selected but --slurm_partition is unset. Set it to your cluster's partition,\n" +
                "  e.g. --slurm_partition global  (or set slurm_partition: <queue> in your params file).\n" +
                "  Find valid partitions with: sinfo -s"
    }
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
