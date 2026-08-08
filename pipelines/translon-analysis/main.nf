#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { validateParameters } from 'plugin/nf-schema'
include { TRANSLON_ANALYSIS } from './workflows/translon_analysis.nf'

workflow {
    validateParameters()
    TRANSLON_ANALYSIS()
}

