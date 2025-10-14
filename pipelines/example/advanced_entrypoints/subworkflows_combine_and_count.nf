include { COMBINE_OUTPUTS } from '../modules/combine_outputs'
include { COUNT_LINES } from '../modules/count_lines'

workflow COMBINE_AND_COUNT {
    take:
    results_channel  // channel: [ val(meta), path(output_file) ]

    main:
    // Group all files by meta.id
    grouped_results = results_channel
        .groupTuple()

    // Combine all outputs into one file
    COMBINE_OUTPUTS(grouped_results)

    // Count lines in the combined file
    COUNT_LINES(COMBINE_OUTPUTS.out.combined)

    // Combine versions
    all_versions = COMBINE_OUTPUTS.out.versions
        .concat(COUNT_LINES.out.versions)

    emit:
    combined = COMBINE_OUTPUTS.out.combined  // channel: [ val(meta), path(combined_file) ]
    counts   = COUNT_LINES.out.counts        // channel: [ val(meta), path(count_file) ]
    versions = all_versions                  // channel: [ path(versions.yml) ]
}
