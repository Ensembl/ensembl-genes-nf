#!/usr/bin/env nextflow

/*
 * Minimal Example Workflow
 *
 * This demonstrates basic subworkflow usage without entry point complexity.
 * It simply runs two subworkflows in sequence.
 */

nextflow.enable.dsl = 2

// Import subworkflows
include { RUN_TOOLS } from '../advanced_entrypoints/subworkflows_run_tools'
include { COMBINE_AND_COUNT } from '../advanced_entrypoints/subworkflows_combine_and_count'

// Parameters
params.outdir = 'results'

workflow {
    // Create input channel
    input_ch = Channel.of(
        [id: 'sample1'],
        [id: 'sample2']
    )
    .map { meta ->
        def dummy_file = file("${workflow.workDir}/dummy_input.txt")
        dummy_file.text = "dummy input"
        [meta, dummy_file]
    }

    // Run the four tools
    RUN_TOOLS(input_ch)

    // Combine and count the results
    COMBINE_AND_COUNT(RUN_TOOLS.out.results)
}
