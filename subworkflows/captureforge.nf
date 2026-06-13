// ── subworkflows/captureforge.nf ───────────────────────────────────────────
// Full Stage 1–10 + visualization pipeline wiring.
include { PARSE_GENE_LIST; DEFINE_TARGETS } from '../modules/target_definition.nf'
include { APPLY_MASKING }                   from '../modules/masking.nf'
include { TILE_CANDIDATES }                 from '../modules/tiling.nf'
include { SPECIFICITY_FILTER }              from '../modules/specificity.nf'
include { COMPOSITION_FILTER }              from '../modules/composition.nf'
include { GAP_RETILE }                      from '../modules/gap_retile.nf'
include { BOOST_DEDUP }                     from '../modules/boost_dedup.nf'
include { SCAFFOLD_SYNTHESIS }              from '../modules/synthesis.nf'
include { AGGREGATE_QC }                    from '../modules/aggregate_qc.nf'
include { EXPORT_PACKAGE }                  from '../modules/export.nf'
include { MAKE_FIGURES; MAKE_DASHBOARD }    from '../modules/visualize.nf'
include { FAIDX }                           from '../modules/faidx.nf'

// optional-input resolver: real file, or a per-slot NO_<slot> sentinel (distinct
// names avoid Nextflow input filename collisions when several are absent).
def optfile(p, slot) {
    p ? file(p, checkIfExists: true) : file("${projectDir}/assets/NO_${slot}")
}

workflow CAPTUREFORGE {
    take:
    gene_list
    gtf
    genome

    main:
    ch_dbsnp   = optfile(params.dbsnp_vcf,        'dbsnp')
    ch_reg     = optfile(params.regulatory_gff,   'regulatory')
    ch_rsid    = optfile(params.rsid_bed,         'rsid')
    ch_rmsk    = optfile(params.rmsk_bed,         'rmsk')
    ch_segdups = optfile(params.segdups_bed,      'segdups')
    ch_simple  = optfile(params.simplerepeat_bed, 'simplerepeat')
    ch_refprov = optfile(params.refs_provenance ?: null, 'refsprov')

    FAIDX(genome)
    fai = FAIDX.out.fai

    // Stage 1
    PARSE_GENE_LIST(gene_list, gtf)
    DEFINE_TARGETS(PARSE_GENE_LIST.out.table, gtf, fai, ch_dbsnp, ch_reg, ch_rsid)
    // Stage 2–3
    APPLY_MASKING(DEFINE_TARGETS.out.raw, genome, fai, ch_rmsk, ch_segdups, ch_simple)
    TILE_CANDIDATES(APPLY_MASKING.out.masked, genome, fai)
    // Stage 4–5
    SPECIFICITY_FILTER(TILE_CANDIDATES.out.fasta, TILE_CANDIDATES.out.bed, fai)
    COMPOSITION_FILTER(SPECIFICITY_FILTER.out.pass, TILE_CANDIDATES.out.fasta)
    // Stage 6–7
    GAP_RETILE(APPLY_MASKING.out.masked, COMPOSITION_FILTER.out.pass,
               COMPOSITION_FILTER.out.fasta, genome, fai)
    BOOST_DEDUP(GAP_RETILE.out.pass, GAP_RETILE.out.fasta)
    // Stage 8
    SCAFFOLD_SYNTHESIS(BOOST_DEDUP.out.tsv)
    // Stage 9 — collect every stage summary for aggregation
    ch_summaries = PARSE_GENE_LIST.out.summary
        .mix(DEFINE_TARGETS.out.summary, APPLY_MASKING.out.summary, TILE_CANDIDATES.out.summary,
             SPECIFICITY_FILTER.out.summary, COMPOSITION_FILTER.out.summary, GAP_RETILE.out.summary,
             BOOST_DEDUP.out.summary, SCAFFOLD_SYNTHESIS.out.summary).collect()
    AGGREGATE_QC(DEFINE_TARGETS.out.raw, APPLY_MASKING.out.masked, BOOST_DEDUP.out.bed, ch_summaries)
    // Stage 10
    EXPORT_PACKAGE(BOOST_DEDUP.out.tsv, BOOST_DEDUP.out.bed, SCAFFOLD_SYNTHESIS.out.tsv,
                   SPECIFICITY_FILTER.out.verdicts, AGGREGATE_QC.out.covered,
                   AGGREGATE_QC.out.metrics, ch_refprov)
    // Visualization — stage every file make_figures reads into one work dir
    ch_vizfiles = AGGREGATE_QC.out.metrics
        .mix(AGGREGATE_QC.out.by_gene, AGGREGATE_QC.out.by_class, COMPOSITION_FILTER.out.verdicts,
             BOOST_DEDUP.out.bed, DEFINE_TARGETS.out.raw, APPLY_MASKING.out.masked,
             SPECIFICITY_FILTER.out.summary, GAP_RETILE.out.summary, BOOST_DEDUP.out.summary).collect()
    MAKE_FIGURES(ch_vizfiles)
    MAKE_DASHBOARD(AGGREGATE_QC.out.metrics, MAKE_FIGURES.out.dir)

    emit:
    covered   = AGGREGATE_QC.out.covered
    baits_csv = EXPORT_PACKAGE.out.csv
    metrics   = AGGREGATE_QC.out.metrics
    dashboard = MAKE_DASHBOARD.out.html
}
