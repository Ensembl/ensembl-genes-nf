/*
 * TRACKHUB_GENERATION SUBWORKFLOW
 * Generate UCSC-compatible track hub from BigWig outputs
 */

include { GENERATE_TRACKHUB } from '../modules/generate_trackhub.nf'

workflow TRACKHUB_GENERATION {
    take:
    bigwigs           // channel: [ meta, bigwig_files ] - per-sample BigWigs
    merged_bigwigs    // channel: [ meta, merged_bigwig_files ] - merged BigWigs
    hub_name          // val: Hub name
    genome            // val: Genome assembly (e.g., 'hg38')
    email             // val: Contact email
    sample_regex      // val: Optional regex for sample ID extraction
    annotation_regex  // val: Optional regex for annotation type extraction

    main:
    // Collect all BigWig files into a single list
    // First, extract just the BigWig paths from the tuples
    all_bigwigs = bigwigs
        .map { _meta, bw_files ->
            bw_files
        }
        .flatten()
        .mix(
            merged_bigwigs
                .map { _meta, bw_files ->
                    bw_files
                }
                .flatten()
        )
        .collect()

    // Generate the track hub
    GENERATE_TRACKHUB(
        all_bigwigs,
        hub_name,
        genome,
        email,
        sample_regex,
        annotation_regex
    )

    emit:
    GENERATE_TRACKHUB.out.trackhub              // path: trackhub directory
}
