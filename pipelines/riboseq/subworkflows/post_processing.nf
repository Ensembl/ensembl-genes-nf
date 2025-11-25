/*
 * POST-PROCESSING SUBWORKFLOW
 * Processes genome BAMs to generate BEDgraph and BigWig files
 * Uses offsets from RiboMetric for P-site correction
 */

include { FILTER_BAM } from '../modules/filter_bam.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FILTERED } from '../modules/samtools_index.nf'
include { BAM_TO_BED } from '../modules/bam_to_bed.nf'
include { MERGE_BEDGRAPHS } from '../modules/merge_bedgraphs.nf'
include { BEDGRAPH_TO_BIGWIG } from '../modules/bedgraph_to_bigwig.nf'
include { BEDGRAPH_TO_BIGWIG as BEDGRAPH_TO_BIGWIG_MERGED } from '../modules/bedgraph_to_bigwig.nf'

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

    // Convert individual BEDgraphs to BigWig
    BEDGRAPH_TO_BIGWIG(
        BAM_TO_BED.out.bedgraph,
        chrom_sizes
    )

    // Collect all bedgraphs per sample for merging
    // Group by meta.id to collect all 4 filtered BAM bedgraphs per sample
    bedgraphs_grouped = BAM_TO_BED.out.bedgraph
        .map { meta, bedgraph ->
            // Flatten the bedgraph list if it's stranded (contains forward/reverse)
            def bg_list = bedgraph instanceof List ? bedgraph : [bedgraph]
            return [meta, bg_list]
        }
        .groupTuple()
        .map { meta, bedgraph_lists ->
            // Flatten all bedgraph lists into a single list per sample
            def all_bedgraphs = bedgraph_lists.flatten()
            return [meta, all_bedgraphs]
        }

    // Merge all bedgraphs per sample (sum coverage across all 4 filtered BAM types)
    MERGE_BEDGRAPHS(bedgraphs_grouped)

    // Convert merged BEDgraphs to BigWig
    BEDGRAPH_TO_BIGWIG_MERGED(
        MERGE_BEDGRAPHS.out.bedgraph,
        chrom_sizes
    )

    emit:
    filtered_bams = FILTER_BAM.out.all_bams                    // tuple: [ meta, bams[] ]
    bedgraphs = BAM_TO_BED.out.bedgraph                        // tuple: [ meta, bedgraph ]
    bigwigs = BEDGRAPH_TO_BIGWIG.out.bigwig                    // tuple: [ meta, bigwig ]
    merged_bedgraphs = MERGE_BEDGRAPHS.out.bedgraph            // tuple: [ meta, merged_bedgraph ]
    merged_bigwigs = BEDGRAPH_TO_BIGWIG_MERGED.out.bigwig      // tuple: [ meta, merged_bigwig ]
}
