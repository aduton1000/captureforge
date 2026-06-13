// ── modules/tiling.nf  (Stage 3) ───────────────────────────────────────────
// Tile candidate baits per-class, then extract their sequences.
process TILE_CANDIDATES {
    tag   "${params.panel_name}"
    label 'medium'
    conda "${projectDir}/env/captureforge.yml"
    publishDir "${params.outdir}/stage3_candidates", mode: params.publish_mode

    input:
    path masked_bed
    path fasta
    path fai

    output:
    path "candidates.bed",      emit: bed
    path "candidates.fasta",    emit: fasta
    path "stage3_summary.json", emit: summary

    script:
    // serialize the per-class tiling policy to JSON for the tiler (no imports under strict syntax)
    def tiling_json = '{' + params.tiling.collect { k, v ->
        "\"${k}\":{\"step_frac\":${v.step_frac},\"even\":${v.even}}"
    }.join(',') + '}'
    """
    tile_candidates.py \\
        --masked ${masked_bed} --fasta ${fasta} --fai ${fai} \\
        --bait-len ${params.bait_len} \\
        --tiling-json '${tiling_json}' \\
        -o candidates.bed --fasta-out candidates.fasta \\
        --summary stage3_summary.json
    """
}
