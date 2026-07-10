/*
 * POST-PROCESSING SUBWORKFLOW
 * Processes genome BAMs to generate BEDgraph and BigWig files
 * Uses offsets from RiboMetric for P-site correction
 * Also creates unique read index and count matrix from collapsed FASTA files
 */

include { FILTER_BAM } from '../modules/filter_bam.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FILTERED } from '../modules/samtools_index.nf'
include { BAM_TO_BED } from '../modules/bam_to_bed.nf'
include { MERGE_BEDGRAPHS } from '../modules/merge_bedgraphs.nf'
include { BEDGRAPH_TO_BIGWIG } from '../modules/bedgraph_to_bigwig.nf'
include { BEDGRAPH_TO_BIGWIG as BEDGRAPH_TO_BIGWIG_MERGED } from '../modules/bedgraph_to_bigwig.nf'
include { MERGE_UNIQUE_READS } from '../modules/merge_unique_reads.nf'

workflow POST_PROCESSING {
    take:
    genome_bam        // tuple: [ meta, bam, bai ]
    offsets           // tuple: [ meta, offsets_file ] - from RiboMetric/RiboWaltz
    chrom_sizes       // path: chromosome sizes file
    collapsed_fastas  // tuple: [ meta, collapsed_fasta ] - for unique read tracking

    main:
    // Filter genome BAM by mapping quality and multiplicity
    FILTER_BAM(
        genome_bam.map { meta, bam, bai -> tuple(meta, bam) }
    )

    // Flatten all 4 filtered BAM outputs into individual emissions with type annotation
    // FILTER_BAM.out.all_bams emits: tuple [ meta, [bam1, bam2, bam3, bam4] ]
    // We need: 4 separate emissions of tuple [ meta_with_type, bam ]
    all_filtered_bams = FILTER_BAM.out.all_bams
        .flatMap { meta, bams ->
            bams.collect { bam ->
                // Extract BAM type from filename (e.g., "SRR123.unique_no_junction.bam" -> "unique_no_junction")
                def bam_name = bam.name
                def bam_type = bam_name.replaceAll(/.*\.([^.]+)\.bam$/, '$1')

                // Add bam_type to meta for tracking
                def meta_with_type = meta + [bam_type: bam_type]
                [meta_with_type, bam]
            }
        }

    // Index all filtered BAMs
    SAMTOOLS_INDEX_FILTERED(all_filtered_bams)
    all_bams_indexed = SAMTOOLS_INDEX_FILTERED.out.bam_and_bai

    // Generate BEDgraph files using offsets
    // Combine each BAM with its corresponding offset file
    // Match by sample ID (meta.id) regardless of bam_type
    bam_with_offsets = all_bams_indexed
        .combine(offsets)
        .filter { bam_meta, bam, bai, offset_meta, offset ->
            bam_meta.id == offset_meta.id
        }
        .map { bam_meta, bam, bai, offset_meta, offset ->
            def tier = offset_meta.track_tier ?: 'selected'
            [bam_meta + [track_tier: tier], bam, bai, offset]
        }

    BAM_TO_BED(bam_with_offsets)

    // Convert individual BEDgraphs to BigWig
    BEDGRAPH_TO_BIGWIG(
        BAM_TO_BED.out.bedgraph,
        chrom_sizes
    )

    // Group bedgraphs by BAM type for separate merged tracks
    // BAM_TO_BED.out.bedgraph emits: tuple [ meta_with_type, [forward.bg, reverse.bg] ]
    bedgraphs_by_type = BAM_TO_BED.out.bedgraph
        .flatMap { meta, bedgraph ->
            // Flatten stranded bedgraphs, keeping bam_type from meta
            def bg_list = bedgraph instanceof List ? bedgraph : [bedgraph]
            bg_list.collect { bg ->
                [meta.track_tier ?: 'selected', meta.bam_type, bg]
            }
        }
        .groupTuple(by: [0, 1])  // Group by track tier and bam_type
        .map { track_tier, bam_type, bedgraph_list ->
            // Create meta with bam_type as id
            def meta = [id: "merged_${track_tier}_${bam_type}", track_tier: track_tier, bam_type: bam_type]
            [meta, bedgraph_list]
        }

    // Merge bedgraphs separately for each BAM type
    MERGE_BEDGRAPHS(bedgraphs_by_type)

    // Convert merged BEDgraphs to BigWig
    BEDGRAPH_TO_BIGWIG_MERGED(
        MERGE_BEDGRAPHS.out.bedgraph,
        chrom_sizes
    )

    // Optional: Collect all collapsed FASTA files and create unique read index
    // Only run if explicitly enabled via params.enable_unique_reads_tracking
    if (params.enable_unique_reads_tracking) {
        all_collapsed_fastas = collapsed_fastas
            .map { meta, fasta -> fasta }
            .collect()

        // Determine mode and previous index
        def unique_reads_mode = params.unique_reads_mode ?: 'per_run'
        def previous_index_dir = params.unique_reads_previous_index ?: null

        // For progressive mode, check if previous index exists
        def previous_index_files = []
        if (unique_reads_mode == 'progressive' && previous_index_dir) {
            def prev_dir = file(previous_index_dir)
            if (prev_dir.exists()) {
                def prev_fasta = file("${previous_index_dir}/unique_reads.fa")
                def prev_matrix = file("${previous_index_dir}/count_matrix.parquet")
                def prev_mapping = file("${previous_index_dir}/read_mapping.json")

                if (prev_fasta.exists() && prev_matrix.exists() && prev_mapping.exists()) {
                    previous_index_files = [prev_fasta, prev_matrix, prev_mapping]
                }
            }
        }

        MERGE_UNIQUE_READS(
            all_collapsed_fastas,
            previous_index_files,
            unique_reads_mode
        )

        // Emit outputs when tracking is enabled
        unique_reads_fasta_out = MERGE_UNIQUE_READS.out.unique_reads_fasta
        count_matrix_out = MERGE_UNIQUE_READS.out.count_matrix
        read_mapping_out = MERGE_UNIQUE_READS.out.read_mapping
        unique_reads_summary_out = MERGE_UNIQUE_READS.out.summary
    } else {
        // Emit empty channels when tracking is disabled
        unique_reads_fasta_out = Channel.empty()
        count_matrix_out = Channel.empty()
        read_mapping_out = Channel.empty()
        unique_reads_summary_out = Channel.empty()
    }

    emit:
    filtered_bams = FILTER_BAM.out.all_bams                    // tuple: [ meta, bams[] ]
    bedgraphs = BAM_TO_BED.out.bedgraph                        // tuple: [ meta, bedgraph ]
    bigwigs = BEDGRAPH_TO_BIGWIG.out.bigwig                    // tuple: [ meta, bigwig ]
    merged_bedgraphs = MERGE_BEDGRAPHS.out.bedgraph            // tuple: [ meta, merged_bedgraph ]
    merged_bigwigs = BEDGRAPH_TO_BIGWIG_MERGED.out.bigwig      // tuple: [ meta, merged_bigwig ]
    unique_reads_fasta = unique_reads_fasta_out                // path: unique_reads.fa (empty if disabled)
    count_matrix = count_matrix_out                            // path: count_matrix.parquet (empty if disabled)
    read_mapping = read_mapping_out                            // path: read_mapping.json (empty if disabled)
    unique_reads_summary = unique_reads_summary_out            // path: processing_summary.json (empty if disabled)
}
