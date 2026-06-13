// ── modules/target_definition.nf  (Stage 1) ────────────────────────────────
// Normalize the gene list, then build the labelled, design-classed target BED.
process PARSE_GENE_LIST {
    tag   "${params.panel_name}"
    label 'small'
    conda "${projectDir}/env/captureforge.yml"
    publishDir "${params.outdir}/stage1_targets", mode: params.publish_mode

    input:
    path gene_list
    path gtf

    output:
    path "targets.tsv",        emit: table
    path "parse_summary.json", emit: summary

    script:
    """
    parse_gene_list.py \\
        --input ${gene_list} \\
        --gtf ${gtf} \\
        --excel-sheet "${params.excel_sheet}" \\
        --default-classes ${params.default_classes.join(',')} \\
        --anchor-flank ${params.anchor_flank} \\
        -o targets.tsv --summary parse_summary.json
    """
}

process DEFINE_TARGETS {
    tag   "${params.panel_name}"
    label 'small'
    conda "${projectDir}/env/captureforge.yml"
    publishDir "${params.outdir}/stage1_targets", mode: params.publish_mode

    input:
    path targets_tsv
    path gtf
    path fai
    path dbsnp          // optional ('NO_FILE' sentinel when absent)
    path regulatory     // optional
    path rsid_bed       // optional pre-resolved rsID BED

    output:
    path "targets.raw.bed",  emit: raw
    path "stage1_summary.json", emit: summary

    script:
    def dbsnp_arg = dbsnp.name != 'NO_FILE'      ? "--dbsnp ${dbsnp}"           : ''
    def reg_arg   = regulatory.name != 'NO_FILE' ? "--regulatory ${regulatory}" : ''
    def rsid_arg  = rsid_bed.name != 'NO_FILE'   ? "--rsid-bed ${rsid_bed}"     : ''
    """
    define_targets.py \\
        --targets ${targets_tsv} --gtf ${gtf} --fai ${fai} \\
        --canonical ${params.canonical_source} \\
        --splice-flank ${params.splice_flank} \\
        --promoter-up ${params.promoter_upstream} --promoter-down ${params.promoter_downstream} \\
        --anchor-flank ${params.anchor_flank} \\
        ${dbsnp_arg} ${reg_arg} ${rsid_arg} \\
        -o targets.raw.bed --summary stage1_summary.json
    """
}
