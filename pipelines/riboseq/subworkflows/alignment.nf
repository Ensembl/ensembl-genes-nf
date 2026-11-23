/*
 * ALIGNMENT SUBWORKFLOW
 * Aligns collapsed reads to genome using STAR
 * STAR generates both genome and transcriptome BAMs
 */

include { STAR_ALIGN } from '../modules/star_align.nf'
include { STAR_ALIGN_GETRPF_HOLDS } from '../modules/star_align_getrpf_holds.nf'
include { SAMTOOLS_SORT } from '../modules/samtools_sort.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_GENOME } from '../modules/samtools_index.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TRANSCRIPTOME } from '../modules/samtools_index.nf'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_GETRPF_HOLDS } from '../modules/samtools_sort.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_GETRPF_HOLDS_GENOME } from '../modules/samtools_index.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_GETRPF_HOLDS_TRANSCRIPTOME } from '../modules/samtools_index.nf'
include { PUBLISH_STAR_ALIGNMENTS } from '../modules/publish_star_alignments.nf'
include { PUBLISH_STAR_ALIGNMENTS as PUBLISH_STAR_ALIGNMENTS_GETRPF_HOLDS } from '../modules/publish_star_alignments.nf'

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

    // Publish only the final indexed STAR alignments using the conventional
    // names consumed by downstream ORF-calling workflows.
    PUBLISH_STAR_ALIGNMENTS(
        SAMTOOLS_INDEX_GENOME.out.bam_and_bai.map { meta, bam, bai -> tuple(meta, bam, bai, 'genome') }
            .mix(
                SAMTOOLS_INDEX_TRANSCRIPTOME.out.bam_and_bai.map { meta, bam, bai -> tuple(meta, bam, bai, 'transcriptome') }
            )
    )


    emit:
    genome_bam = SAMTOOLS_INDEX_GENOME.out.bam_and_bai           // tuple: [ meta, bam, bai ]
    transcriptome_bam = SAMTOOLS_INDEX_TRANSCRIPTOME.out.bam_and_bai                // tuple: [ meta, bam, bai ]
    logs = STAR_ALIGN.out.log                                    // tuple: [ meta, log ]
}

/*
 * Diagnostic alignment of raw reads withheld by getRPF. The Local-STAR BAMs
 * are published separately and are never connected to analysis or QC gating.
 */
workflow ALIGNMENT_GETRPF_HOLDS {
    take:
    raw_reads          // tuple: [ meta, raw_fastq.gz ]
    star_index         // path: STAR index directory
    gtf                // path: GTF annotation file

    main:
    diagnostic_reads = raw_reads.map { meta, reads ->
        tuple(meta + [alignment_purpose: 'getrpf_hold_diagnostic'], reads)
    }

    STAR_ALIGN_GETRPF_HOLDS(
        diagnostic_reads,
        star_index,
        gtf
    )

    SAMTOOLS_INDEX_GETRPF_HOLDS_GENOME(STAR_ALIGN_GETRPF_HOLDS.out.bam)

    SAMTOOLS_SORT_GETRPF_HOLDS(STAR_ALIGN_GETRPF_HOLDS.out.transcriptome_bam)
    SAMTOOLS_INDEX_GETRPF_HOLDS_TRANSCRIPTOME(SAMTOOLS_SORT_GETRPF_HOLDS.out.bam)

    PUBLISH_STAR_ALIGNMENTS_GETRPF_HOLDS(
        SAMTOOLS_INDEX_GETRPF_HOLDS_GENOME.out.bam_and_bai.map { meta, bam, bai -> tuple(meta, bam, bai, 'genome') }
            .mix(
                SAMTOOLS_INDEX_GETRPF_HOLDS_TRANSCRIPTOME.out.bam_and_bai.map { meta, bam, bai -> tuple(meta, bam, bai, 'transcriptome') }
            )
    )

    emit:
    genome_bam = SAMTOOLS_INDEX_GETRPF_HOLDS_GENOME.out.bam_and_bai
    transcriptome_bam = SAMTOOLS_INDEX_GETRPF_HOLDS_TRANSCRIPTOME.out.bam_and_bai
    logs = STAR_ALIGN_GETRPF_HOLDS.out.log
}
