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
    // Combine each BAM with its corresponding offset file
    // Since all 4 BAMs per sample have the same meta.id, we use combine + filter
    bam_with_offsets = all_bams_indexed
        .combine(offsets)
        .filter { bam_meta, bam, bai, offset_meta, offset ->
            bam_meta.id == offset_meta.id
        }
        .map { bam_meta, bam, bai, offset_meta, offset ->
            [bam_meta, bam, bai, offset]
        }

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

    // Collect all collapsed FASTA files and create unique read index
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

    emit:
    filtered_bams = FILTER_BAM.out.all_bams                    // tuple: [ meta, bams[] ]
    bedgraphs = BAM_TO_BED.out.bedgraph                        // tuple: [ meta, bedgraph ]
    bigwigs = BEDGRAPH_TO_BIGWIG.out.bigwig                    // tuple: [ meta, bigwig ]
    merged_bedgraphs = MERGE_BEDGRAPHS.out.bedgraph            // tuple: [ meta, merged_bedgraph ]
    merged_bigwigs = BEDGRAPH_TO_BIGWIG_MERGED.out.bigwig      // tuple: [ meta, merged_bigwig ]
    unique_reads_fasta = MERGE_UNIQUE_READS.out.unique_reads_fasta     // path: unique_reads.fa
    count_matrix = MERGE_UNIQUE_READS.out.count_matrix                 // path: count_matrix.parquet
    read_mapping = MERGE_UNIQUE_READS.out.read_mapping                 // path: read_mapping.json
    unique_reads_summary = MERGE_UNIQUE_READS.out.summary              // path: processing_summary.json
}
