#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SUBWORKFLOW_EXAMPLE } from './workflows/subworkflow_example.nf'

// Parameters
params.outdir = 'results'

workflow {
    SUBWORKFLOW_EXAMPLE()
}
