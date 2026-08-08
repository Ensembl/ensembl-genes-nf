/*
 * ANALYSIS SUBWORKFLOW
 * Runs RiboMetric for QC and offset calculation when an annotation is available.
 * RiboWaltz is kept as a fallback, or as an explicit offset source when
 * params.ribometric_use_ribowaltz_offsets = true.
 */

include { RIBOMETRIC } from '../modules/ribometric.nf'
include { RIBOWALTZ } from '../modules/ribowaltz.nf'

workflow ANALYSIS {
    take:
    transcriptome_bam      // tuple: [ meta, bam, bai ]
    ribometric_annotation  // path: RiboMetric annotation file (optional)
    gtf                    // path: GTF annotation file
    fasta                  // path: Reference genome FASTA (for RiboWaltz fallback)

    main:
    def use_ribometric = params.ribometric_annotation || params.run_organism_setup || params.auto_build_indices

    if (use_ribometric) {
        if (params.ribometric_use_ribowaltz_offsets) {
            // Optional legacy mode: calculate offsets with RiboWaltz, then pass
            // them into RiboMetric.
            gtf_ch = gtf.map { value -> [[ id: 'reference' ], value] }
            fasta_ch = fasta.map { value -> [[ id: 'reference' ], value] }

            RIBOWALTZ(
                transcriptome_bam,
                gtf_ch,
                fasta_ch
            )

            // Join RiboWaltz offsets with transcriptome BAM on sample ID
            // RiboWaltz best_offset: tuple [ meta, offset_file ]
            // transcriptome_bam: tuple [ meta, bam, bai ]
            ribometric_input = transcriptome_bam
                .join(RIBOWALTZ.out.best_offset)
                .map { meta, bam, bai, offset_file ->
                    tuple(meta, bam, bai, offset_file)
                }

            RIBOMETRIC(
                ribometric_input.map { meta, bam, bai, _offset -> tuple(meta, bam, bai) },
                ribometric_annotation,
                ribometric_input.map { _meta, _bam, _bai, offset -> offset }
            )

            psite_offsets_ch = RIBOWALTZ.out.psite_offsets
            best_offset_ch = RIBOWALTZ.out.best_offset
            psite_table_ch = RIBOWALTZ.out.psite_table
            cds_coverage_ch = RIBOWALTZ.out.cds_coverage
            codon_rpf_ch = RIBOWALTZ.out.codon_rpf
            codon_psite_ch = RIBOWALTZ.out.codon_psite
            offset_plots_ch = RIBOWALTZ.out.offset_plots
            qc_plots_ch = RIBOWALTZ.out.qc_plots
        } else {
            // Default path: use RiboMetric's internal calculation and avoid
            // the more resource-intensive RiboWaltz task.
            RIBOMETRIC(
                transcriptome_bam,
                ribometric_annotation,
                file('NO_OFFSET_FILE')  // Placeholder for optional input
            )

            psite_offsets_ch = channel.empty()
            best_offset_ch = channel.empty()
            psite_table_ch = channel.empty()
            cds_coverage_ch = channel.empty()
            codon_rpf_ch = channel.empty()
            codon_psite_ch = channel.empty()
            offset_plots_ch = channel.empty()
            qc_plots_ch = channel.empty()
        }

        ribometric_html_ch = RIBOMETRIC.out.html
        ribometric_json_ch = RIBOMETRIC.out.json
        ribometric_csv_ch = RIBOMETRIC.out.csv
        ribometric_offsets_audit_ch = RIBOMETRIC.out.offsets_audit
        offsets_ch = RIBOMETRIC.out.offsets
    } else {
        // No RiboMetric annotation: fall back to RiboWaltz offsets so
        // downstream BEDgraph generation can still proceed.
        gtf_ch = gtf.map { value -> [[ id: 'reference' ], value] }
        fasta_ch = fasta.map { value -> [[ id: 'reference' ], value] }

        RIBOWALTZ(
            transcriptome_bam,
            gtf_ch,
            fasta_ch
        )

        ribometric_html_ch = channel.empty()
        ribometric_json_ch = channel.empty()
        ribometric_csv_ch = channel.empty()
        ribometric_offsets_audit_ch = channel.empty()
        offsets_ch = RIBOWALTZ.out.best_offset

        psite_offsets_ch = RIBOWALTZ.out.psite_offsets
        best_offset_ch = RIBOWALTZ.out.best_offset
        psite_table_ch = RIBOWALTZ.out.psite_table
        cds_coverage_ch = RIBOWALTZ.out.cds_coverage
        codon_rpf_ch = RIBOWALTZ.out.codon_rpf
        codon_psite_ch = RIBOWALTZ.out.codon_psite
        offset_plots_ch = RIBOWALTZ.out.offset_plots
        qc_plots_ch = RIBOWALTZ.out.qc_plots
    }

    emit:
    // RiboMetric outputs
    ribometric_html = ribometric_html_ch
    ribometric_json = ribometric_json_ch
    ribometric_csv = ribometric_csv_ch
    ribometric_offsets_audit = ribometric_offsets_audit_ch

    // RiboWaltz outputs
    psite_offsets = psite_offsets_ch                 // tuple: [ meta, tsv.gz ]
    best_offset = best_offset_ch                     // tuple: [ meta, txt ]
    psite_table = psite_table_ch                     // tuple: [ meta, tsv.gz ]
    cds_coverage = cds_coverage_ch                   // tuple: [ meta, tsv.gz ]
    codon_rpf = codon_rpf_ch                         // tuple: [ meta, tsv.gz ]
    codon_psite = codon_psite_ch                     // tuple: [ meta, tsv.gz ]
    offset_plots = offset_plots_ch                   // tuple: [ meta, pdfs ]
    qc_plots = qc_plots_ch                           // tuple: [ meta, pdfs ]

    offsets = offsets_ch
}
