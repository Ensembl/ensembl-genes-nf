#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { validateParameters } from 'plugin/nf-schema'

// Validate parameters against schema
validateParameters()

include { SUBWORKFLOW_EXAMPLE } from './workflows/subworkflow_example.nf'

workflow {
    SUBWORKFLOW_EXAMPLE()
}