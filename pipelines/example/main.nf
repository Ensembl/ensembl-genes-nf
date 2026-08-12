#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

nextflow.enable.strict = true

include { validateParameters } from 'plugin/nf-schema'

include { SUBWORKFLOW_EXAMPLE } from './workflows/subworkflow_example.nf'

workflow {
    validateParameters()
    SUBWORKFLOW_EXAMPLE(file(params.input))
}
