// ── modules/synthesis.nf  (Stage 8) ────────────────────────────────────────
process SCAFFOLD_SYNTHESIS {
    tag   "${params.panel_name}"
    label 'small'
    conda "${projectDir}/env/captureforge.yml"
    publishDir "${params.outdir}/stage8_synthesis", mode: params.publish_mode

    input:
    path baits_tsv

    output:
    path "synthesis.fasta",     emit: fasta
    path "synthesis.tsv",       emit: tsv
    path "stage8_summary.json", emit: summary

    script:
    """
    scaffold_synthesis.py \\
        --baits-tsv ${baits_tsv} --format ${params.bait_format} \\
        --t7 ${params.t7_promoter} --handle5 ${params.pcr_handle_5} --handle3 ${params.pcr_handle_3} \\
        --max-len ${params.synth_max_len} --strand ${params.bait_strand} \\
        --prefix synthesis --summary stage8_summary.json
    """
}
