/*
 * ANALYSIS SUBWORKFLOW
 * Runs both RiboMetric and RiboWaltz for complementary QC and offset calculation
 * - RiboMetric: Fast QC and offset calculation
 * - RiboWaltz: Detailed P-site analysis and comprehensive profiles
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
    // Run RiboMetric if annotation provided
    if (ribometric_annotation) {
        RIBOMETRIC(
            transcriptome_bam,
            ribometric_annotation
        )
    }

    // Prepare GTF and FASTA channels that can be reused for all samples
    gtf_ch = Channel.value([[ id: 'reference' ], gtf])
    fasta_ch = Channel.value([[ id: 'reference' ], fasta])

    // Run RiboWaltz on transcriptome BAM
    RIBOWALTZ(
        transcriptome_bam,
        gtf_ch,
        fasta_ch
    )

    emit:
    // RiboMetric outputs
    ribometric_html = ribometric_annotation ? RIBOMETRIC.out.html : Channel.empty()
    ribometric_json = ribometric_annotation ? RIBOMETRIC.out.json : Channel.empty()
    ribometric_csv = ribometric_annotation ? RIBOMETRIC.out.csv : Channel.empty()
    ribometric_offsets = ribometric_annotation ? RIBOMETRIC.out.offsets : Channel.empty()

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
    offsets = RIBOWALTZ.out.best_offset
}
