// ── modules/faidx.nf ────────────────────────────────────────────────────────
// Ensure a .fai exists for the (possibly large) genome FASTA. Idempotent.
process FAIDX {
    tag   "${fasta.name}"
    label 'small'
    conda "${projectDir}/env/captureforge.yml"

    input:
    path fasta

    output:
    path "${fasta.name}.fai", emit: fai

    script:
    """
    if [ -f "${fasta}.fai" ]; then
        cp -P "${fasta}.fai" "${fasta.name}.fai" 2>/dev/null || cp "${fasta}.fai" "${fasta.name}.fai"
    else
        samtools faidx ${fasta}
    fi
    """
}
