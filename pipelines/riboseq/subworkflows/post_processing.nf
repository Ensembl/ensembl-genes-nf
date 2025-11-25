/*
 * POST-PROCESSING SUBWORKFLOW
 * Processes genome BAMs to generate BEDgraph and BigWig files
 * Uses offsets from RiboMetric for P-site correction
 */

include { FILTER_BAM } from '../modules/filter_bam.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FILTERED } from '../modules/samtools_index.nf'
include { BAM_TO_BED } from '../modules/bam_to_bed.nf'
include { BEDGRAPH_TO_BIGWIG } from '../modules/bedgraph_to_bigwig.nf'

workflow POST_PROCESSING {
    take:
    genome_bam        // tuple: [ meta, bam, bai ]
    offsets           // tuple: [ meta, offsets_file ] - from RiboMetric/RiboWaltz
    chrom_sizes       // path: chromosome sizes file

    main:
    // Filter genome BAM by mapping quality and multiplicity
    FILTER_BAM(
        genome_bam.map { meta, bam, bai -> tuple(meta, bam) }
    )

    // Collect all 4 filtered BAM outputs into a single channel
    // Each emits: tuple [ meta, bam ]
    all_filtered_bams = FILTER_BAM.out.unique_no_junction
        .concat(FILTER_BAM.out.unique_with_junction)
        .concat(FILTER_BAM.out.multi_no_junction)
        .concat(FILTER_BAM.out.multi_with_junction)

    // Index all filtered BAMs
    SAMTOOLS_INDEX_FILTERED(all_filtered_bams)
    all_bams_indexed = SAMTOOLS_INDEX_FILTERED.out.bam_and_bai

    // Generate BEDgraph files using offsets
    // Join filtered BAMs with offsets
    bam_with_offsets = all_bams_indexed
        .join(offsets)

    BAM_TO_BED(bam_with_offsets)

    // Convert BEDgraph to BigWig
    BEDGRAPH_TO_BIGWIG(
        BAM_TO_BED.out.bedgraph,
        chrom_sizes
    )

    emit:
    filtered_bams = FILTER_BAM.out.all_bams           // tuple: [ meta, bams[] ]
    bedgraphs = BAM_TO_BED.out.bedgraph               // tuple: [ meta, bedgraph ]
    bigwigs = BEDGRAPH_TO_BIGWIG.out.bigwig           // tuple: [ meta, bigwig ]
}
