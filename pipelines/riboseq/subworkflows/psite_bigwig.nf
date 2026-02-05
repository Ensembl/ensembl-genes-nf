/*
 * PSITE_BIGWIG SUBWORKFLOW
 * Takes paired genome and transcriptome BAMs and produces P-site bigwig files
 *
 * Flow:
 *   1. RIBOWALTZ: Calculate P-site offsets from transcriptome BAM
 *   2. FILTER_BAM: Filter genome BAM by mapping quality/multiplicity (optional)
 *   3. BAM_TO_BED: Convert genome BAM to bedgraph using P-site offsets
 *   4. BEDGRAPH_TO_BIGWIG: Convert bedgraph to bigwig
 */

include { RIBOWALTZ } from '../modules/ribowaltz.nf'
include { FILTER_BAM } from '../modules/filter_bam.nf'
include { SAMTOOLS_INDEX } from '../modules/samtools_index.nf'
include { BAM_TO_BED } from '../modules/bam_to_bed.nf'
include { BEDGRAPH_TO_BIGWIG } from '../modules/bedgraph_to_bigwig.nf'

workflow PSITE_BIGWIG {
    take:
    genome_bam         // tuple: [ meta, bam, bai ]
    transcriptome_bam  // tuple: [ meta, bam, bai ]
    gtf                // path: GTF annotation file
    fasta              // path: Reference genome FASTA
    chrom_sizes        // path: Chromosome sizes file

    main:
    // Prepare GTF and FASTA channels with metadata for RiboWaltz
    gtf_ch = gtf.map { [[ id: 'reference' ], it] }
    fasta_ch = fasta.map { [[ id: 'reference' ], it] }

    // Step 1: Calculate P-site offsets from transcriptome BAM using RiboWaltz
    RIBOWALTZ(
        transcriptome_bam,
        gtf_ch,
        fasta_ch
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

    // Step 3: Combine BAMs with their corresponding P-site offsets
    // Match by sample ID (meta.id) regardless of bam_type
    bam_with_offsets = bams_to_process
        .combine(RIBOWALTZ.out.best_offset)
        .filter { bam_meta, bam, bai, offset_meta, offset ->
            bam_meta.id == offset_meta.id
        }
        .map { bam_meta, bam, bai, offset_meta, offset ->
            [bam_meta, bam, bai, offset]
        }

    // Step 4: Convert BAMs to BEDgraph using P-site offsets
    BAM_TO_BED(bam_with_offsets)

    // Step 5: Convert BEDgraphs to BigWig
    BEDGRAPH_TO_BIGWIG(
        BAM_TO_BED.out.bedgraph,
        chrom_sizes
    )

    emit:
    // RiboWaltz outputs
    psite_offsets = RIBOWALTZ.out.psite_offsets     // tuple: [ meta, tsv.gz ]
    best_offset = RIBOWALTZ.out.best_offset          // tuple: [ meta, txt ]
    psite_table = RIBOWALTZ.out.psite_table          // tuple: [ meta, tsv.gz ]
    qc_plots = RIBOWALTZ.out.qc_plots                // tuple: [ meta, pdfs ]

    // Bedgraph and BigWig outputs
    bedgraphs = BAM_TO_BED.out.bedgraph              // tuple: [ meta, bedgraph ]
    bigwigs = BEDGRAPH_TO_BIGWIG.out.bigwig          // tuple: [ meta, bigwig ]
}
