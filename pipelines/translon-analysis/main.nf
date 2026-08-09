#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params {
    tools: String = 'all-wave1'
    samplesheet: String = null
    min_caller_agreement: Integer = 2
    run_characterisation: Boolean = true
    threads_ribocode: Integer = 4
    threads_ribotricer: Integer = 4
    threads_ribotaper: Integer = 4
    threads_orfquant: Integer = 4
    threads_rpbp: Integer = 4
}

include { validateParameters } from 'plugin/nf-schema'
include { TRANSLON_ANALYSIS } from './workflows/translon_analysis.nf'

workflow {
    validateParameters()
    TRANSLON_ANALYSIS()
}
