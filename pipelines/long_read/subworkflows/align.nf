/*
 * ALIGN SUBWORKFLOW
 * minimap2 alignment + samtools sort/index for one long-read sample
 */

include { MINIMAP2_ALIGN } from '../modules/minimap2_align.nf'
include { SAMTOOLS_SORT  } from '../modules/samtools_sort.nf'
include { SAMTOOLS_INDEX } from '../modules/samtools_index.nf'

workflow ALIGN {
    take:
    reads_ch   // [ meta, fastq ]
    index_ch   // [ meta, mmi  ]

    main:
    ch_versions = Channel.empty()

    MINIMAP2_ALIGN(
        reads_ch,
        index_ch,
        true,    // bam_format
        'bai',   // bam_index_extension
        false,   // cigar_paf_format
        true     // cigar_bam — handles CIGAR >65535 ops for long reads
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    SAMTOOLS_SORT(
        MINIMAP2_ALIGN.out.bam,
        [[],[],[]],  // no reference needed for BAM output
        'bai'
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    // Join sorted BAM with its index for downstream consumption
    bam_with_index = SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_INDEX.out.index)

    emit:
    bam      = SAMTOOLS_SORT.out.bam    // [ meta, bam ]
    bai      = SAMTOOLS_INDEX.out.index // [ meta, bai ]
    bam_bai  = bam_with_index           // [ meta, bam, bai ]
    versions = ch_versions
}
