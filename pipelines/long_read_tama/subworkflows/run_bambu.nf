nextflow.enable.dsl = 2

include { BAMBU_DISCOVERY } from '../modules/bambu.nf'
include { GTF_TO_BED12 } from '../modules/gtf_to_bed12.nf'

workflow RUN_BAMBU {
    take:
    aligned_bams
    reference

    main:
    analysis = aligned_bams
        .map { meta, bam, bai -> tuple(meta, [bam], [bai], 'accession', params.isoquant_data_type) }
    if (params.bambu_scope in ['cohort', 'both']) {
        cohort = aligned_bams
            .map { _meta, bam, bai -> tuple('cohort', bam, bai) }
            .groupTuple()
            .map { _key, bams, bais -> tuple([id: params.cohort_id, scope: 'cohort'], bams, bais, 'cohort', params.isoquant_data_type) }
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
    versions = BAMBU_DISCOVERY.out.versions.mix(GTF_TO_BED12.out.versions)
}
