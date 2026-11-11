nextflow.enable.dsl = 2

include { ACQUIRE_READS } from '../modules/acquire_reads.nf'
include { TRIM_CAGE_READS } from '../modules/trim_cage_reads.nf'
include { ALIGN_CAGE_READS } from '../modules/align_cage_reads.nf'

workflow PROCESS_READS {
    take:
    samples
    indexes

    main:
    ch_acquired = ACQUIRE_READS(samples).reads
    ch_trimmed = TRIM_CAGE_READS(ch_acquired).reads
    ch_trimmed_keyed = ch_trimmed.map { meta, reads -> tuple(meta.assembly, meta, reads) }
    ch_indexes_keyed = indexes.map { meta, index, excluded -> tuple(meta.assembly, index, excluded) }
    ch_joined = ch_trimmed_keyed.join(ch_indexes_keyed)
        .map { _assembly, meta, reads, index, excluded -> tuple(meta, reads, index, excluded) }
    ch_bam = ALIGN_CAGE_READS(ch_joined).bam

    emit:
    bam = ch_bam
}
