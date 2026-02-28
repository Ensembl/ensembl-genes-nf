/*
 * ANALYSIS SUBWORKFLOW
 * Runs both RiboMetric and RiboWaltz for complementary QC and offset calculation
 * - RiboMetric: Fast QC and offset calculation
 * - RiboWaltz: Detailed P-site analysis and comprehensive profiles
 *
 * Supports passing RiboWaltz offsets to RiboMetric when params.ribometric_use_ribowaltz_offsets = true
 */

include { RIBOMETRIC } from '../modules/ribometric.nf'
include { RIBOWALTZ } from '../modules/ribowaltz.nf'

workflow ANALYSIS {
    take:
    transcriptome_bam      // tuple: [ meta, bam, bai ]
    ribometric_annotation  // path: RiboMetric annotation file (optional)
    gtf                    // path: GTF annotation file
    fasta                  // path: Reference genome FASTA (for RiboWaltz)

    main:
    // Prepare GTF and FASTA channels with metadata for RiboWaltz
    // gtf and fasta are already channels, just need to add metadata
    gtf_ch = gtf.map { [[ id: 'reference' ], it] }
    fasta_ch = fasta.map { [[ id: 'reference' ], it] }

    // Run RiboWaltz on transcriptome BAM (always runs first when passing offsets)
    RIBOWALTZ(
        transcriptome_bam,
        gtf_ch,
        fasta_ch
    )

    // Run RiboMetric if annotation provided
    if (ribometric_annotation) {
        // Determine offset input for RiboMetric
        if (params.ribometric_use_ribowaltz_offsets) {
            // Join RiboWaltz offsets with transcriptome BAM on sample ID
            // RiboWaltz best_offset: tuple [ meta, offset_file ]
            // transcriptome_bam: tuple [ meta, bam, bai ]
            ribometric_input = transcriptome_bam
                .join(RIBOWALTZ.out.best_offset)
                .map { meta, bam, bai, offset_file ->
                    [ meta, bam, bai, offset_file ]
                }

            RIBOMETRIC(
                ribometric_input.map { meta, bam, bai, offset -> [ meta, bam, bai ] },
                ribometric_annotation,
                ribometric_input.map { meta, bam, bai, offset -> offset }
            )
        } else {
            // No external offset file - use RiboMetric's internal calculation
            RIBOMETRIC(
                transcriptome_bam,
                ribometric_annotation,
                file('NO_OFFSET_FILE')  // Placeholder for optional input
            )
        }
    }

    emit:
    // RiboMetric outputs
    ribometric_html = ribometric_annotation ? RIBOMETRIC.out.html : Channel.empty()
    ribometric_json = ribometric_annotation ? RIBOMETRIC.out.json : Channel.empty()
    ribometric_csv = ribometric_annotation ? RIBOMETRIC.out.csv : Channel.empty()

    // RiboWaltz outputs
    psite_offsets = RIBOWALTZ.out.psite_offsets     // tuple: [ meta, tsv.gz ]
    best_offset = RIBOWALTZ.out.best_offset          // tuple: [ meta, txt ]
    psite_table = RIBOWALTZ.out.psite_table          // tuple: [ meta, tsv.gz ]
    cds_coverage = RIBOWALTZ.out.cds_coverage        // tuple: [ meta, tsv.gz ]
    codon_rpf = RIBOWALTZ.out.codon_rpf              // tuple: [ meta, tsv.gz ]
    codon_psite = RIBOWALTZ.out.codon_psite          // tuple: [ meta, tsv.gz ]
    offset_plots = RIBOWALTZ.out.offset_plots        // tuple: [ meta, pdfs ]
    qc_plots = RIBOWALTZ.out.qc_plots                // tuple: [ meta, pdfs ]

    // Use RiboWaltz best_offset for downstream processing
    offsets = RIBOMETRIC.out.offsets
}
