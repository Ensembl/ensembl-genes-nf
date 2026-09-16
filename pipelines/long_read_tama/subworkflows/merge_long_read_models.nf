nextflow.enable.dsl = 2

include { TAMA_MERGE } from '../modules/tama_merge.nf'
include { TMERGE } from '../modules/tmerge.nf'

workflow MERGE_LONG_READ_MODELS {
    take:
    beds

    main:
    // beds: tuple val(cohort_id), path(list of TAMA BED files)
    if (params.merge_tool == 'tmerge') {
        TMERGE(beds)
        bed = TMERGE.out.bed
        filelist = TMERGE.out.filelist
        reports = TMERGE.out.gene_report.mix(TMERGE.out.merge_report).mix(TMERGE.out.trans_report)
        versions = TMERGE.out.versions
    } else {
        TAMA_MERGE(beds)
        bed = TAMA_MERGE.out.bed
        filelist = TAMA_MERGE.out.filelist
        reports = TAMA_MERGE.out.gene_report.mix(TAMA_MERGE.out.merge_report).mix(TAMA_MERGE.out.trans_report)
        versions = TAMA_MERGE.out.versions
    }

    emit:
    bed
    filelist
    reports
    versions
}
