nextflow.enable.dsl = 2

include { EXTRACT_CTSS } from '../modules/extract_ctss.nf'

workflow BUILD_CTSS {
    take:
    bam

    main:
    ch_ctss = EXTRACT_CTSS(bam, file("${projectDir}/bin/extract_ctss.sh")).ctss

    emit:
    ctss = ch_ctss
}
