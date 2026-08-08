#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { validateParameters } from 'plugin/nf-schema'
include { TRANSLON_CHARACTERISATION } from './workflows/translon_characterisation.nf'

workflow {
    validateParameters()

    TRANSLON_CHARACTERISATION(
        file(params.intervals, checkIfExists: true),
        file(params.gencode_gff3, checkIfExists: true),
        file(params.translation_verdicts, checkIfExists: true),
        file(params.proteome_fasta, checkIfExists: true),
        file(params.genome_fasta, checkIfExists: true)
    )
}
