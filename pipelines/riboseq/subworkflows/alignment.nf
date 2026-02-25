/*
 * ALIGNMENT SUBWORKFLOW
 * Aligns collapsed reads to genome using STAR
 * STAR generates both genome and transcriptome BAMs
 */

include { STAR_ALIGN } from '../modules/star_align.nf'
include { SAMTOOLS_SORT } from '../modules/samtools_sort.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_GENOME } from '../modules/samtools_index.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TRANSCRIPTOME } from '../modules/samtools_index.nf'

workflow ALIGNMENT {
    take:
    collapsed_reads   // tuple: [ meta, collapsed_fasta ]
    star_index        // path: STAR genome index directory
    gtf               // path: GTF annotation file

    main:
    // Run STAR alignment
    // STAR outputs genome BAM (sorted) and transcriptome BAM (unsorted)
    STAR_ALIGN(
        collapsed_reads,
        star_index,
        gtf
    )

    // Index genome BAM (already sorted by coordinate)
    SAMTOOLS_INDEX_GENOME(STAR_ALIGN.out.bam)

    // Sort transcriptome BAM before indexing
    SAMTOOLS_SORT(STAR_ALIGN.out.transcriptome_bam)

    // Index sorted transcriptome BAM
    SAMTOOLS_INDEX_TRANSCRIPTOME(SAMTOOLS_SORT.out.bam)


    emit:
    genome_bam = SAMTOOLS_INDEX_GENOME.out.bam_and_bai           // tuple: [ meta, bam, bai ]
    transcriptome_bam = SAMTOOLS_INDEX_TRANSCRIPTOME.out.bam_and_bai                // tuple: [ meta, bam, bai ]
    logs = STAR_ALIGN.out.log                                    // tuple: [ meta, log ]
}
