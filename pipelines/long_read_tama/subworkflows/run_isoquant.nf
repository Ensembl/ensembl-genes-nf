nextflow.enable.dsl = 2

include { ISOQUANT_ANNOTATION_FREE; ISOQUANT_REFERENCE_GUIDED } from '../modules/isoquant.nf'
include { NORMALISE_ISOQUANT_OUTPUTS } from '../modules/normalise_isoquant_outputs.nf'
include { AUDIT_ISOQUANT_OUTPUTS } from '../modules/audit_isoquant_outputs.nf'
include { ISOQUANT_TO_BED12 } from '../modules/isoquant_to_bed12.nf'
include { REPORT_UNSUPPORTED_MODEL_INPUTS as REPORT_UNSUPPORTED_ISOQUANT_INPUTS } from '../modules/report_unsupported_model_inputs.nf'

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
    // IsoQuant has no raw-subread data_type. Run PacBio raw-subread BAMs with
    // its PacBio algorithm settings, while preserving raw_subread provenance
    // in source_read_metadata and never changing the biological read label.
    isoquant_type_for = { meta ->
        meta.isoquant_data_type ?: (meta.technology == 'pacbio' && meta.read_representation == 'raw_subread' ? 'pacbio_ccs' : null)
    }
    eligible_bams = aligned_bams.filter { meta, _bam, _bai ->
        isoquant_type_for(meta) in ['pacbio_ccs', 'nanopore']
    }
    unsupported_bams = aligned_bams.filter { meta, _bam, _bai ->
        !(isoquant_type_for(meta) in ['pacbio_ccs', 'nanopore'])
    }
    REPORT_UNSUPPORTED_ISOQUANT_INPUTS(unsupported_bams.map { meta, _bam, _bai -> tuple(meta, 'isoquant') })
    common = eligible_bams.map { meta, bam, bai ->
        def read_metadata = "${meta.id}:${meta.technology}:${meta.read_representation}:${meta.preprocessing_state}:${meta.observation_unit}"
        def data_type = isoquant_type_for(meta)
        tuple(meta + [source_read_metadata: read_metadata, isoquant_data_type: data_type], [bam], [bai], 'accession', data_type, analysis, large_output, read_group, params.isoquant_args ?: '', params.isoquant_report_novel_unspliced, params.isoquant_check_canonical)
    }
    if (params.isoquant_scope in ['cohort', 'both']) {
        cohort = eligible_bams
            .map { meta, bam, bai -> tuple(isoquant_type_for(meta), meta, bam, bai) }
            .groupTuple()
            .map { data_type, metas, bams, bais ->
                def source_metadata = metas.sort { left, right -> left.id <=> right.id }.collect { source_meta -> "${source_meta.id}:${source_meta.technology}:${source_meta.read_representation}:${source_meta.preprocessing_state}:${source_meta.observation_unit}" }.join(';')
                tuple([id: "${params.cohort_id}_${data_type}", scope: 'cohort', isoquant_data_type: data_type, source_accessions: metas.collect { source_meta -> source_meta.id }.sort().join(','), source_read_metadata: source_metadata], bams, bais, 'cohort', data_type, analysis, large_output, read_group, params.isoquant_args ?: '', params.isoquant_report_novel_unspliced, params.isoquant_check_canonical)
            }
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
    reports = AUDIT_ISOQUANT_OUTPUTS.out.manifest.mix(AUDIT_ISOQUANT_OUTPUTS.out.status).mix(REPORT_UNSUPPORTED_ISOQUANT_INPUTS.out.report)
    versions = versions.mix(NORMALISE_ISOQUANT_OUTPUTS.out.versions)
        .mix(AUDIT_ISOQUANT_OUTPUTS.out.versions)
        .mix(ISOQUANT_TO_BED12.out.versions)
        .mix(REPORT_UNSUPPORTED_ISOQUANT_INPUTS.out.versions)
}
