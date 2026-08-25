#!/usr/bin/env nextflow

include { validateParameters; samplesheetToList } from 'plugin/nf-schema'
include { IMPORT_REFSEQ } from './workflows/import_refseq.nf'

workflow {

    validateParameters()

    sample_ch = Channel
        .fromList(samplesheetToList(params.input_csv, "${projectDir}/assets/schema_input.json"))
        .map { gcf, species ->
            [
                id     : gcf,
                species: species
            ]
        }

    def settings = new groovy.json.JsonSlurper().parse(file(params.server_settings))
    if (!settings.db_read_user)
        throw new IllegalArgumentException("server_settings must contain db_read_user for taxonomy loading")

    Channel
        .value(
            tuple(
                settings.db_host,
                settings.db_port,
                settings.db_user,
                settings.db_password,
                settings.db_read_user
            )
        )
        .set { db_config_ch }

    IMPORT_REFSEQ(
        sample_ch,
        db_config_ch
    )
}
