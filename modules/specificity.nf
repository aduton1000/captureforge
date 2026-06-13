// ── modules/specificity.nf  (Stage 4) ──────────────────────────────────────
// Whole-genome off-target detection. The genome index / BLAST db is referenced
// in place via params.specificity_db (absolute prefix) — NOT staged, since it
// is multi-GB. This is the heavy stage (label 'index'); on HPC it gets a big slot.
process SPECIFICITY_FILTER {
    tag   "${params.panel_name}"
    label 'index'
    conda "${projectDir}/env/captureforge.yml"
    publishDir "${params.outdir}/stage4_specificity", mode: params.publish_mode

    input:
    path candidates_fasta
    path candidates_bed

    output:
    path "candidates.specpass.bed", emit: pass
    path "candidates.spec.tsv",     emit: verdicts
    path "stage4_summary.json",     emit: summary

    script:
    if (!params.specificity_db)
        error "SPECIFICITY_FILTER requires --specificity_db (aligner index or BLAST db prefix)"
    """
    specificity_filter.py \\
        --candidates ${candidates_fasta} --candidates-bed ${candidates_bed} \\
        --aligner ${params.specificity_aligner} --db ${params.specificity_db} \\
        --min-identity ${params.offtarget_min_identity} \\
        --min-lenfrac ${params.offtarget_min_lenfrac} \\
        --max-offtarget ${params.offtarget_max_hits} \\
        --threads ${task.cpus} \\
        -o candidates.specpass.bed --verdicts candidates.spec.tsv \\
        --summary stage4_summary.json
    """
}
