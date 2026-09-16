nextflow.enable.dsl = 2

include { BUILD_MINIMAP2_INDEX } from '../modules/build_minimap_index.nf'
include { MINIMAP2_ALIGN } from '../modules/minimap2_align.nf'
include { SAMTOOLS_SORT_INDEX } from '../modules/samtools_sort_index.nf'
include { ALIGNMENT_QC } from '../modules/alignment_qc.nf'

workflow ALIGN_LONG_READS {
    take:
    reads
    reference

    main:
    // reads: tuple val(meta), path(reads.fastq.gz)
    // bam: tuple val(meta), path(sorted.bam), path(sorted.bam.bai)
    build_versions = channel.empty()
    if (params.minimap_index) {
        minimap_index = channel.value(file(params.minimap_index, checkIfExists: true))
    } else {
        BUILD_MINIMAP2_INDEX(reference)
        // The reference is a single value, so the single generated index is
        // also reusable by every alignment task when paired with the read stream.
        minimap_index = BUILD_MINIMAP2_INDEX.out.index
        build_versions = BUILD_MINIMAP2_INDEX.out.versions
    }
    MINIMAP2_ALIGN(reads, minimap_index)
    SAMTOOLS_SORT_INDEX(MINIMAP2_ALIGN.out.sam)
    ALIGNMENT_QC(SAMTOOLS_SORT_INDEX.out.bam)

    emit:
    bam = SAMTOOLS_SORT_INDEX.out.bam
    stats = ALIGNMENT_QC.out.stats
    versions = build_versions.mix(MINIMAP2_ALIGN.out.versions).mix(SAMTOOLS_SORT_INDEX.out.versions).mix(ALIGNMENT_QC.out.versions)
}
