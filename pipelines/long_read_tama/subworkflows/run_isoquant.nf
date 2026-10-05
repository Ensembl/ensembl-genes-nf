nextflow.enable.dsl = 2

include { ISOQUANT_ANNOTATION_FREE; ISOQUANT_REFERENCE_GUIDED } from '../modules/isoquant.nf'
include { NORMALISE_ISOQUANT_OUTPUTS } from '../modules/normalise_isoquant_outputs.nf'
include { AUDIT_ISOQUANT_OUTPUTS } from '../modules/audit_isoquant_outputs.nf'
include { ISOQUANT_TO_BED12 } from '../modules/isoquant_to_bed12.nf'

workflow RUN_ISOQUANT {
    take:
    aligned_bams
    reference

    main:
    // aligned_bams: tuple val(meta), path(sorted.bam), path(sorted.bam.bai)
    // bed: tuple val(meta), val(scope), path(isoquant_models.bed), cohort only
    analysis = params.isoquant_analysis instanceof List ? params.isoquant_analysis : [params.isoquant_analysis]
    large_output = params.isoquant_large_output instanceof List ? params.isoquant_large_output : [params.isoquant_large_output]
    read_group = params.isoquant_read_group instanceof List ? params.isoquant_read_group : [params.isoquant_read_group]
    common = aligned_bams.map { meta, bam, bai -> tuple(meta, [bam], [bai], 'accession', params.isoquant_data_type, analysis, large_output, read_group, params.isoquant_args ?: '', params.isoquant_report_novel_unspliced, params.isoquant_check_canonical) }
    if (params.isoquant_scope in ['cohort', 'both']) {
        cohort = aligned_bams
            .map { _meta, bam, bai -> tuple('cohort', bam, bai) }
            .groupTuple()
            .map { _key, bams, bais -> tuple([id: params.cohort_id, scope: 'cohort'], bams, bais, 'cohort', params.isoquant_data_type, analysis, large_output, read_group, params.isoquant_args ?: '', params.isoquant_report_novel_unspliced, params.isoquant_check_canonical) }
    } else {
        cohort = channel.empty()
    }

    jobs = params.isoquant_scope == 'accession' ? common : params.isoquant_scope == 'cohort' ? cohort : common.mix(cohort)
    if (params.isoquant_mode == 'annotation_free') {
        ISOQUANT_ANNOTATION_FREE(jobs, reference)
        isoquant_native = ISOQUANT_ANNOTATION_FREE.out.isoquant_native
        versions = ISOQUANT_ANNOTATION_FREE.out.versions
    } else {
        genedb = channel.value(file(params.isoquant_genedb, checkIfExists: true))
        ISOQUANT_REFERENCE_GUIDED(jobs, reference, genedb)
        isoquant_native = ISOQUANT_REFERENCE_GUIDED.out.isoquant_native
        versions = ISOQUANT_REFERENCE_GUIDED.out.versions
    }
    NORMALISE_ISOQUANT_OUTPUTS(isoquant_native)
    AUDIT_ISOQUANT_OUTPUTS(NORMALISE_ISOQUANT_OUTPUTS.out.products)
    ISOQUANT_TO_BED12(NORMALISE_ISOQUANT_OUTPUTS.out.gtf)

    emit:
    products = NORMALISE_ISOQUANT_OUTPUTS.out.products
    gtf = NORMALISE_ISOQUANT_OUTPUTS.out.gtf
    bed = ISOQUANT_TO_BED12.out.bed.filter { _meta, scope, _bed -> scope == 'cohort' }
        .map { meta, scope, bed -> tuple(meta, 'isoquant', scope, bed) }
    accession_bed = ISOQUANT_TO_BED12.out.bed.filter { _meta, scope, _bed -> scope == 'accession' }
        .map { meta, scope, bed -> tuple(meta, 'isoquant', scope, bed) }
    reports = AUDIT_ISOQUANT_OUTPUTS.out.manifest.mix(AUDIT_ISOQUANT_OUTPUTS.out.status)
    versions = versions.mix(NORMALISE_ISOQUANT_OUTPUTS.out.versions)
        .mix(AUDIT_ISOQUANT_OUTPUTS.out.versions)
        .mix(ISOQUANT_TO_BED12.out.versions)
}
