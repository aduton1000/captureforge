// ── modules/gap_retile.nf  (Stage 6) ───────────────────────────────────────
// Gap analysis + adaptive re-tiling. Re-uses Stage 4/5 tools internally, so it
// also needs the specificity DB (label 'index' = heavy, like Stage 4).
process GAP_RETILE {
    tag   "${params.panel_name}"
    label 'index'
    conda "${projectDir}/env/captureforge.yml"
    publishDir "${params.outdir}/stage6_retile", mode: params.publish_mode

    input:
    path masked_bed
    path pass_bed
    path pass_fasta
    path fasta
    path fai

    output:
    path "retiled.pass.bed",   emit: pass
    path "retiled.pass.fasta", emit: fasta
    path "residual_gaps.bed",  emit: residual
    path "stage6_summary.json", emit: summary

    script:
    if (!params.specificity_db)
        error "GAP_RETILE requires --specificity_db"
    """
    gap_retile.py \\
        --masked ${masked_bed} --passing-bed ${pass_bed} --passing-fasta ${pass_fasta} \\
        --fasta ${fasta} --fai ${fai} --db ${params.specificity_db} --aligner ${params.specificity_aligner} \\
        --bait-len ${params.bait_len} --max-passes ${params.retile_max_passes} \\
        --coverage-target ${params.coverage_target} --threads ${task.cpus} \\
        --min-identity ${params.offtarget_min_identity} --min-lenfrac ${params.offtarget_min_lenfrac} \\
        --max-offtarget ${params.offtarget_max_hits} \\
        --gc-min ${params.gc_min} --gc-max ${params.gc_max} \\
        --tm-min ${params.tm_min} --tm-max ${params.tm_max} --homopolymer-max ${params.homopolymer_max} \\
        -o retiled.pass.bed --fasta-out retiled.pass.fasta --residual residual_gaps.bed \\
        --summary stage6_summary.json
    """
}
