nextflow.enable.dsl = 2

include { TAMA_MERGE } from '../modules/tama_merge.nf'

workflow MERGE_LONG_READ_MODELS {
    take:
    beds
    filelist_builder

    main:
    // beds: tuple val(cohort_id), path(list of TAMA BED files)
    TAMA_MERGE(beds, filelist_builder)

    emit:
    bed = TAMA_MERGE.out.bed
    filelist = TAMA_MERGE.out.filelist
    reports = TAMA_MERGE.out.gene_report.mix(TAMA_MERGE.out.merge_report).mix(TAMA_MERGE.out.trans_report)
    versions = TAMA_MERGE.out.versions
}
