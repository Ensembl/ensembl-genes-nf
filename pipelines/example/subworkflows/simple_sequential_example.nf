/*
 * Simple Sequential Subworkflow Example
 *
 * This shows a subworkflow where processes run sequentially
 * (output of one becomes input to the next)
 */

include { COMBINE_OUTPUTS } from '../modules/combine_outputs'
include { COUNT_LINES } from '../modules/count_lines'

workflow SIMPLE_SEQUENTIAL_EXAMPLE {
    take:
    input_files  // channel: [ val(meta), path(files) ]

    main:
    // First step: combine files
    COMBINE_OUTPUTS(input_files.groupTuple())

    // Second step: count lines in combined file
    COUNT_LINES(COMBINE_OUTPUTS.out.combined)

    emit:
    combined = COMBINE_OUTPUTS.out.combined  // channel: [ val(meta), path(combined_file) ]
    counts   = COUNT_LINES.out.counts        // channel: [ val(meta), path(count_file) ]
    versions = COMBINE_OUTPUTS.out.versions.concat(COUNT_LINES.out.versions)
}
