// ── modules/export.nf  (Stage 10) ──────────────────────────────────────────
process EXPORT_PACKAGE {
    tag   "${params.panel_name}"
    label 'small'
    conda "${projectDir}/env/captureforge.yml"
    publishDir "${params.outdir}/order_package", mode: params.publish_mode

    input:
    path baits_tsv
    path baits_bed
    path synthesis_tsv
    path spec_tsv
    path covered_bed
    path metrics
    path refs_provenance      // 'NO_FILE' sentinel allowed

    output:
    path "baits.fasta",      emit: fasta
    path "baits.csv",        emit: csv
    path "order_generic.tsv", emit: order
    path "MANIFEST",         emit: manifest
    path "provenance.json",  emit: provenance

    script:
    def refs_arg = refs_provenance.name.startsWith('NO_') ? '' : "--refs-provenance ${refs_provenance}"
    """
    export_package.py \\
        --baits-tsv ${baits_tsv} --baits-bed ${baits_bed} --synthesis-tsv ${synthesis_tsv} \\
        --spec-tsv ${spec_tsv} --format ${params.bait_format} --outdir . \\
        --genome "${params.genome_fasta}" --specificity-db "${params.specificity_db}" \\
        --manifest-files ${covered_bed},${metrics} ${refs_arg} \\
        --summary stage10_summary.json
    """
}
