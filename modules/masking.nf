// ── modules/masking.nf  (Stage 2) ──────────────────────────────────────────
// Per-class masking policy + maskability QC.
process APPLY_MASKING {
    tag   "${params.panel_name}"
    label 'medium'
    conda "${projectDir}/env/captureforge.yml"
    publishDir "${params.outdir}/stage2_masking", mode: params.publish_mode

    input:
    path raw_bed
    path fasta
    path fai
    path rmsk           // optional ('NO_FILE' sentinel)
    path segdups        // optional
    path simplerepeat   // optional

    output:
    path "targets.masked.bed",        emit: masked
    path "targets.lost_to_mask.bed",  emit: lost
    path "targets.masked.str_core.bed", emit: str_core, optional: true
    path "stage2_summary.json",       emit: summary

    script:
    def rmsk_arg = rmsk.name    != 'NO_FILE' ? "--rmsk ${rmsk}"       : ''
    def seg_arg  = segdups.name != 'NO_FILE' ? "--segdups ${segdups}" : ''
    def sr_arg   = simplerepeat.name != 'NO_FILE' ? "--simplerepeat ${simplerepeat}" : ''
    """
    apply_masking.py \\
        --raw ${raw_bed} --fasta ${fasta} \\
        ${rmsk_arg} ${seg_arg} ${sr_arg} \\
        -o targets.masked.bed --lost targets.lost_to_mask.bed \\
        --summary stage2_summary.json
    """
}
