nextflow.enable.dsl = 2

include { BAMBU_DISCOVERY } from '../modules/bambu.nf'
include { GTF_TO_BED12 } from '../modules/gtf_to_bed12.nf'
include { REPORT_UNSUPPORTED_MODEL_INPUTS as REPORT_UNSUPPORTED_BAMBU_INPUTS } from '../modules/report_unsupported_model_inputs.nf'

workflow RUN_BAMBU {
    take:
    aligned_bams
    reference

    main:
    eligible_bams = aligned_bams.filter { meta, _bam, _bai -> meta.isoquant_data_type in ['pacbio_ccs', 'nanopore'] }
    unsupported_bams = aligned_bams.filter { meta, _bam, _bai -> !(meta.isoquant_data_type in ['pacbio_ccs', 'nanopore']) }
    REPORT_UNSUPPORTED_BAMBU_INPUTS(unsupported_bams.map { meta, _bam, _bai -> tuple(meta, 'bambu') })
    analysis = eligible_bams
        .map { meta, bam, bai ->
            def read_metadata = "${meta.id}:${meta.technology}:${meta.read_representation}:${meta.preprocessing_state}:${meta.observation_unit}"
            tuple(meta + [source_read_metadata: read_metadata], [bam], [bai], 'accession', meta.isoquant_data_type)
        }
    if (params.bambu_scope in ['cohort', 'both']) {
        cohort = eligible_bams
            .map { meta, bam, bai -> tuple(meta.isoquant_data_type, meta, bam, bai) }
            .groupTuple()
            .map { data_type, metas, bams, bais ->
                def source_metadata = metas.sort { left, right -> left.id <=> right.id }.collect { source_meta -> "${source_meta.id}:${source_meta.technology}:${source_meta.read_representation}:${source_meta.preprocessing_state}:${source_meta.observation_unit}" }.join(';')
                tuple([id: "${params.cohort_id}_${data_type}", scope: 'cohort', isoquant_data_type: data_type, source_accessions: metas.collect { source_meta -> source_meta.id }.sort().join(','), source_read_metadata: source_metadata], bams, bais, 'cohort', data_type)
            }
    } else {
        cohort = channel.empty()
    }
    jobs = params.bambu_scope == 'accession' ? analysis : params.bambu_scope == 'cohort' ? cohort : analysis.mix(cohort)
    BAMBU_DISCOVERY(jobs, reference)
    GTF_TO_BED12(BAMBU_DISCOVERY.out.gtf.map { meta, scope, gtf -> tuple(meta, 'bambu', scope, gtf) })

    emit:
    products = BAMBU_DISCOVERY.out.products
    gtf = BAMBU_DISCOVERY.out.gtf
    bed = GTF_TO_BED12.out.bed.filter { _meta, _backend, scope, _bed -> scope == 'cohort' }
        .map { meta, _backend, scope, bed -> tuple(meta, 'bambu', scope, bed) }
    accession_bed = GTF_TO_BED12.out.bed.filter { _meta, _backend, scope, _bed -> scope == 'accession' }
        .map { meta, _backend, scope, bed -> tuple(meta, 'bambu', scope, bed) }
    versions = BAMBU_DISCOVERY.out.versions.mix(GTF_TO_BED12.out.versions).mix(REPORT_UNSUPPORTED_BAMBU_INPUTS.out.versions)
    reports = REPORT_UNSUPPORTED_BAMBU_INPUTS.out.report
}
