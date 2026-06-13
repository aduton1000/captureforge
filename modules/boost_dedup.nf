// ── modules/boost_dedup.nf  (Stage 7) ──────────────────────────────────────
// Boost copies over flagged regions + collapse identical sequences.
process BOOST_DEDUP {
    tag   "${params.panel_name}"
    label 'medium'
    conda "${projectDir}/env/captureforge.yml"
    publishDir "${params.outdir}/stage7_boost_dedup", mode: params.publish_mode

    input:
    path retiled_bed
    path retiled_fasta

    output:
    path "baits.dedup.tsv",   emit: tsv
    path "baits.dedup.bed",   emit: bed
    path "baits.dedup.fasta", emit: fasta
    path "stage7_summary.json", emit: summary

    script:
    """
    boost_dedup.py \\
        --bed ${retiled_bed} --fasta ${retiled_fasta} \\
        --per-oligo-price ${params.per_oligo_price} \\
        --pool-tiers ${params.pool_tiers.join(',')} \\
        --prefix baits.dedup --summary stage7_summary.json
    """
}
