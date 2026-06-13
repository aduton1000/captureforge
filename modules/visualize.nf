// ── modules/visualize.nf  (figures + self-contained HTML dashboard) ─────────
process MAKE_FIGURES {
    tag   "${params.panel_name}"
    label 'medium'
    conda "${projectDir}/env/captureforge.yml"
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    path workfiles   // staged: BEDs, *.tsv, *.json, metrics.json (collected)

    output:
    path "figures/*.png", emit: png
    path "figures/*.svg", emit: svg
    path "figures",       emit: dir

    script:
    """
    mkdir -p figures
    make_figures.py --workdir . --figdir figures --tracks-per-page ${params.tracks_per_page ?: 6}
    """
}

process MAKE_DASHBOARD {
    tag   "${params.panel_name}"
    label 'small'
    conda "${projectDir}/env/captureforge.yml"
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    path metrics
    path figures

    output:
    path "qc_dashboard.html", emit: html

    script:
    """
    make_dashboard.py --metrics ${metrics} --figdir ${figures} \\
        --panel "${params.panel_name}" --bait-format ${params.bait_format} \\
        --generated "\$(date -u +%Y-%m-%dT%H:%M:%SZ)" -o qc_dashboard.html
    """
}
