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

    // Flatten all 4 filtered BAM outputs into individual emissions
    // FILTER_BAM.out.all_bams emits: tuple [ meta, [bam1, bam2, bam3, bam4] ]
    // We need: 4 separate emissions of tuple [ meta, bam ]
    all_filtered_bams = FILTER_BAM.out.all_bams
        .flatMap { meta, bams ->
            bams.collect { bam ->
                [meta, bam]
            }
        }

    // Index all filtered BAMs
    SAMTOOLS_INDEX_FILTERED(all_filtered_bams)
    all_bams_indexed = SAMTOOLS_INDEX_FILTERED.out.bam_and_bai

    // Generate BEDgraph files using offsets
    // Join filtered BAMs with offsets
    // Now each of the 4 BAMs will be joined with the same offset file
    bam_with_offsets = all_bams_indexed
        .join(offsets)

    BAM_TO_BED(bam_with_offsets)

    // Convert individual BEDgraphs to BigWig
    BEDGRAPH_TO_BIGWIG(
        BAM_TO_BED.out.bedgraph,
        chrom_sizes
    )

    // Collect ALL bedgraphs across ALL samples for grand total merge
    // BAM_TO_BED.out.bedgraph emits: tuple [ meta, [forward.bg, reverse.bg] ] per BAM
    // We want to merge everything into a single aggregate track
    all_bedgraphs = BAM_TO_BED.out.bedgraph
        .flatMap { meta, bedgraph ->
            // Flatten stranded bedgraphs into separate files
            def bg_list = bedgraph instanceof List ? bedgraph : [bedgraph]
            bg_list
        }
        .collect()
        .map { bedgraph_list ->
            // Create a single meta for the merged output with id 'all_merged'
            def meta = [id: 'all_merged']
            [meta, bedgraph_list]
        }

    // Merge ALL bedgraphs (sum coverage across all samples and all filtering types)
    MERGE_BEDGRAPHS(all_bedgraphs)

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
