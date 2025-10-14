#!/usr/bin/env nextflow

/*
 * Simple Workflow - Bare Minimum Example
 *
 * Shows the absolute simplest workflow structure:
 * - Create input
 * - Run one process
 * - Done
 */

nextflow.enable.dsl = 2

include { TOOL_A } from '../modules/tool_a'

// Parameters with defaults
params.outdir = 'simple_results'

workflow {
    // Create a simple input channel
    Channel.of(
        [id: 'sample1'],
        [id: 'sample2']
    )
    .map { meta ->
        def input_file = file("${workflow.workDir}/input.txt")
        input_file.text = "input data"
        [meta, input_file]
    }
    .set { input_ch }

    // Run one process
    TOOL_A(input_ch)
}
