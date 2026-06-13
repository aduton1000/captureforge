// ── modules/composition.nf  (Stage 5) ──────────────────────────────────────
// GC / homopolymer / low-complexity / Tm / optional ΔG filtering.
process COMPOSITION_FILTER {
    tag   "${params.panel_name}"
    label 'medium'
    conda "${projectDir}/env/captureforge.yml"
    publishDir "${params.outdir}/stage5_composition", mode: params.publish_mode

    input:
    path specpass_bed
    path candidates_fasta

    output:
    path "composition.pass.bed",   emit: pass
    path "composition.pass.fasta", emit: fasta
    path "composition.tsv",        emit: verdicts
    path "stage5_summary.json",    emit: summary

    script:
    def vienna = params.use_vienna ? "--use-vienna --deltaG-min ${params.deltaG_min}" : ''
    """
    composition_filter.py \\
        --pass-bed ${specpass_bed} --fasta ${candidates_fasta} \\
        --gc-min ${params.gc_min} --gc-max ${params.gc_max} \\
        --gc-max-promoter ${params.gc_max_promoter} \\
        --homopolymer-max ${params.homopolymer_max} \\
        --tm-min ${params.tm_min} --tm-max ${params.tm_max} \\
        ${vienna} \\
        --prefix composition --verdicts composition.tsv --summary stage5_summary.json
    """
}
