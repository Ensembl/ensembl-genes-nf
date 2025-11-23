/*
 * ANALYSIS SUBWORKFLOW
 * Runs RiboMetric (or RiboWaltz) for QC and offset calculation
 */

include { RIBOMETRIC } from '../modules/ribometric.nf'

workflow ANALYSIS {
    take:
    transcriptome_bam   // tuple: [ meta, bam, bai ]
    annotation          // path: RiboMetric annotation file

    main:
    // Run RiboMetric on transcriptome BAM
    RIBOMETRIC(
        transcriptome_bam,
        annotation
    )

    emit:
    html = RIBOMETRIC.out.html        // tuple: [ meta, html ]
    json = RIBOMETRIC.out.json        // tuple: [ meta, json ]
    csv = RIBOMETRIC.out.csv          // tuple: [ meta, csv ]
    // TODO: Extract offsets from RiboMetric output for post-processing
    // For now, this is a placeholder
    offsets = RIBOMETRIC.out.json.map { meta, json ->
        // In the future, parse JSON to extract offsets
        tuple(meta, json)
    }
}
