/*
 * Minimal Subworkflow Example
 *
 * This shows the basic structure of a subworkflow:
 * - Takes input channels
 * - Runs one or more processes
 * - Emits output channels
 */

include { TOOL_A } from '../modules/tool_a'
include { TOOL_B } from '../modules/tool_b'

workflow MINIMAL_SUBWORKFLOW_EXAMPLE {
    take:
    input_channel  // channel: [ val(meta), path(input_file) ]

    main:
    // Run processes
    TOOL_A(input_channel)
    TOOL_B(input_channel)

    // Combine outputs if needed
    all_results = TOOL_A.out.results
        .concat(TOOL_B.out.results)

    emit:
    results  = all_results      // channel: [ val(meta), path(output_file) ]
    versions = TOOL_A.out.versions.concat(TOOL_B.out.versions)
}
