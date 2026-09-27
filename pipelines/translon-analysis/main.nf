#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params {
    tools: String = 'all'
    samplesheet: String = null
    min_caller_agreement: Integer = 2
    run_consensus: Boolean = true
    run_characterisation: Boolean = true
    merge_inputs: Boolean = false
    merge_group: String = 'all'
    skip_fastq_tools: Boolean = false
    ribotricer_external_offsets: Boolean = true
    partition_mode: String = 'transcriptome'
    partition_count: Integer = 1
    partition_window_size: Integer = 5000000
    partition_padding: Integer = 0
    partition_fai: String = null
    translon_analysis_bin: String = "${projectDir}/bin"
}

include { validateParameters } from 'plugin/nf-schema'
include { TRANSLON_ANALYSIS } from './workflows/translon_analysis.nf'

workflow {
    validateParameters()
    TRANSLON_ANALYSIS()
}
