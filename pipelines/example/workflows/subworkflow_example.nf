#!/usr/bin/env nextflow

/*
 * Subworkflow Example
 *
 * Demonstrates chaining two subworkflows together.
 * - First: Run 2 tools in parallel
 * - Second: Combine outputs and count lines
 */

nextflow.enable.dsl = 2

include { MINIMAL_SUBWORKFLOW_EXAMPLE } from '../subworkflows/minimal_subworkflow_example'
include { SIMPLE_SEQUENTIAL_EXAMPLE } from '../subworkflows/simple_sequential_example'

params.outdir = 'results'

workflow SUBWORKFLOW_EXAMPLE {
    // Create input
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

    // First subworkflow: run 2 tools in parallel
    MINIMAL_SUBWORKFLOW_EXAMPLE(input_ch)

    // Second subworkflow: combine and count
    SIMPLE_SEQUENTIAL_EXAMPLE(MINIMAL_SUBWORKFLOW_EXAMPLE.out.results)
}
