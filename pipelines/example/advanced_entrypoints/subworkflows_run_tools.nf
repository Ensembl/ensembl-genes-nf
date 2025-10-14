include { TOOL_A } from '../modules/tool_a'
include { TOOL_B } from '../modules/tool_b'
include { TOOL_C } from '../modules/tool_c'
include { TOOL_D } from '../modules/tool_d'

workflow RUN_TOOLS {
    take:
    input_channel  // channel: [ val(meta), path(input_file) ]

    main:
    // Run all four tools in parallel
    TOOL_A(input_channel)
    TOOL_B(input_channel)
    TOOL_C(input_channel)
    TOOL_D(input_channel)

    // Combine all results into a single channel
    all_results = TOOL_A.out.results
        .concat(TOOL_B.out.results)
        .concat(TOOL_C.out.results)
        .concat(TOOL_D.out.results)

    // Combine all versions
    all_versions = TOOL_A.out.versions
        .concat(TOOL_B.out.versions)
        .concat(TOOL_C.out.versions)
        .concat(TOOL_D.out.versions)

    emit:
    results  = all_results   // channel: [ val(meta), path(output_file) ]
    versions = all_versions  // channel: [ path(versions.yml) ]
}
