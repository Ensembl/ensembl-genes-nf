/*
 * PSITE_BIGWIG SUBWORKFLOW
 * Takes paired genome and transcriptome BAMs and produces P-site bigwig files
 *
 * Flow:
 *   1. RIBOMETRIC: Calculate A-site offsets from transcriptome BAM
 *   2. FILTER_BAM: Filter genome BAM by mapping quality/multiplicity (optional)
 *   3. BAM_TO_BED: Convert genome BAM to bedgraph using offsets
 *   4. BEDGRAPH_TO_BIGWIG: Convert bedgraph to bigwig
 *
 * Note: requires RiboMetric built with --output-offsets support.
 */

include { SAMTOOLS_SORT as SORT_TRANSCRIPTOME }      from '../modules/samtools_sort.nf'
include { SAMTOOLS_INDEX as INDEX_TRANSCRIPTOME }     from '../modules/samtools_index.nf'
include { RIBOMETRIC }         from '../modules/ribometric.nf'
include { FILTER_BAM }         from '../modules/filter_bam.nf'
include { SAMTOOLS_INDEX }     from '../modules/samtools_index.nf'
include { BAM_TO_BED }         from '../modules/bam_to_bed.nf'
include { BEDGRAPH_TO_BIGWIG } from '../modules/bedgraph_to_bigwig.nf'

workflow PSITE_BIGWIG {
    take:
    genome_bam             // tuple: [ meta, bam, bai ]
    transcriptome_bam      // tuple: [ meta, bam, bai ]
    ribometric_annotation  // path: RiboMetric annotation TSV file
    chrom_sizes            // path: Chromosome sizes file

    main:
    SORT_TRANSCRIPTOME(
        transcriptome_bam.map { meta, bam, bai -> tuple(meta, bam) }
    )

    INDEX_TRANSCRIPTOME(SORT_TRANSCRIPTOME.out.bam)
    // Step 1: Calculate A-site offsets from transcriptome BAM using RiboMetric
    RIBOMETRIC(
        INDEX_TRANSCRIPTOME.out.bam_and_bai,
        ribometric_annotation,
        file('NO_OFFSET_FILE')  // Use internal offset calculation
    )

    // Step 2: Optionally filter genome BAMs
    if (params.filter_bams != false) {
        FILTER_BAM(
            genome_bam.map { meta, bam, bai -> tuple(meta, bam) }
        )

        // Flatten filtered BAMs into individual emissions with type annotation
        all_filtered_bams = FILTER_BAM.out.all_bams
            .flatMap { meta, bams ->
                bams.collect { bam ->
                    def bam_name = bam.name
                    def bam_type = bam_name.replaceAll(/.*\.([^.]+)\.bam$/, '$1')
                    def meta_with_type = meta + [bam_type: bam_type]
                    [meta_with_type, bam]
                }
            }

        // Index filtered BAMs
        SAMTOOLS_INDEX(all_filtered_bams)
        bams_to_process = SAMTOOLS_INDEX.out.bam_and_bai
    } else {
        // Use original genome BAMs without filtering
        bams_to_process = genome_bam.map { meta, bam, bai ->
            def meta_with_type = meta + [bam_type: 'all']
            [meta_with_type, bam, bai]
        }
    }

    // Step 3: Combine BAMs with their corresponding offsets
    // Normalize the key to the sample name before '.Aligned' to match across
    // sortedByCoord (genome) and toTranscriptome (transcriptome) naming conventions
    bam_with_offsets = bams_to_process
        .map { meta, bam, bai -> [meta.id.split('\\.Aligned')[0], meta, bam, bai] }
        .join(
            RIBOMETRIC.out.offsets.map { meta, offset -> [meta.id.split('\\.Aligned')[0], offset] }
        )
        .map { key, bam_meta, bam, bai, offset -> [bam_meta, bam, bai, offset] }

    // Step 4: Convert BAMs to BEDgraph using offsets
    BAM_TO_BED(bam_with_offsets)

    // Step 5: Convert BEDgraphs to BigWig
    BEDGRAPH_TO_BIGWIG(
        BAM_TO_BED.out.bedgraph,
        chrom_sizes
    )

    emit:
    // RiboMetric outputs
    ribometric_html = RIBOMETRIC.out.html    // tuple: [ meta, html ]
    ribometric_json = RIBOMETRIC.out.json    // tuple: [ meta, json ]
    ribometric_csv  = RIBOMETRIC.out.csv     // tuple: [ meta, csv ]
    offsets         = RIBOMETRIC.out.offsets // tuple: [ meta, tsv ]

    // Bedgraph and BigWig outputs
    bedgraphs = BAM_TO_BED.out.bedgraph     // tuple: [ meta, bedgraph ]
    bigwigs   = BEDGRAPH_TO_BIGWIG.out.bigwig // tuple: [ meta, bigwig ]
}
