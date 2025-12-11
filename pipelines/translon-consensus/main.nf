#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { validateParameters } from 'plugin/nf-schema'

// Validate parameters against schema
validateParameters()

include { TRANSLON_CONSENSUS } from './workflows/translon_consensus.nf'

workflow {
    TRANSLON_CONSENSUS()
}