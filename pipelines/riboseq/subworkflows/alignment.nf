/*
 * ALIGNMENT SUBWORKFLOW
 * Aligns collapsed reads to genome using STAR
 * STAR generates both genome and transcriptome BAMs
 */

include { STAR_ALIGN } from '../modules/star_align.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_GENOME } from '../modules/samtools_index.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TRANSCRIPTOME } from '../modules/samtools_index.nf'

workflow ALIGNMENT {
    take:
    collapsed_reads   // tuple: [ meta, collapsed_fasta ]
    star_index        // path: STAR genome index directory
    gtf               // path: GTF annotation file

    main:
    // Run STAR alignment
    // STAR outputs both genome BAM and transcriptome BAM (if enabled)
    STAR_ALIGN(
        collapsed_reads,
        star_index,
        gtf
    )

    // Index genome BAM
    SAMTOOLS_INDEX_GENOME(STAR_ALIGN.out.bam)

    // Index transcriptome BAM if it exists
    if (params.save_star_transcriptome_bam) {
        SAMTOOLS_INDEX_TRANSCRIPTOME(STAR_ALIGN.out.transcriptome_bam)
        transcriptome_bam_indexed = SAMTOOLS_INDEX_TRANSCRIPTOME.out.bam_and_bai
    } else {
        transcriptome_bam_indexed = Channel.empty()
    }

    emit:
    genome_bam = SAMTOOLS_INDEX_GENOME.out.bam_and_bai           // tuple: [ meta, bam, bai ]
    transcriptome_bam = transcriptome_bam_indexed                // tuple: [ meta, bam, bai ]
    logs = STAR_ALIGN.out.log                                    // tuple: [ meta, log ]
}
