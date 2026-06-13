// ── modules/aggregate_qc.nf  (Stage 9) ─────────────────────────────────────
// Aggregate metrics, emit the AS-BUILT handoff BED, run QC cross-checks.
process AGGREGATE_QC {
    tag   "${params.panel_name}"
    label 'medium'
    conda "${projectDir}/env/captureforge.yml"
    publishDir "${params.outdir}/stage9_qc", mode: params.publish_mode

    input:
    path raw_bed
    path masked_bed
    path final_bed
    path summaries          // all stageN_summary.json, staged into the work dir

    output:
    path "final_covered_targets.bed", emit: covered
    path "metrics.json",              emit: metrics
    path "coverage_by_gene.tsv",      emit: by_gene
    path "coverage_by_class.tsv",     emit: by_class
    path "stage9_summary.json",       emit: summary

    script:
    """
    mkdir -p summaries && cp -f *summary*.json summaries/ 2>/dev/null || true
    aggregate_qc.py \\
        --raw ${raw_bed} --masked ${masked_bed} --final-bed ${final_bed} \\
        --summaries-dir summaries \\
        --per-oligo-price ${params.per_oligo_price} --pool-tiers ${params.pool_tiers.join(',')} \\
        --covered-vs-raw-min ${params.covered_vs_raw_min} --bait-len ${params.bait_len} \\
        -o final_covered_targets.bed --metrics metrics.json \\
        --by-gene coverage_by_gene.tsv --by-class coverage_by_class.tsv --summary stage9_summary.json
    """
}
