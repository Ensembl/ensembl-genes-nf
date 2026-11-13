#!/usr/bin/env nextflow

include { validateParameters; samplesheetToList } from 'plugin/nf-schema'
include { IMPORT_REFSEQ } from './workflows/import_refseq.nf'

workflow {

    validateParameters()

    sample_ch = channel
        .fromList(samplesheetToList(params.input_csv, "${projectDir}/assets/schema_input.json"))
        .map { gcf, species ->
            [
                id     : gcf,
                species: species
            ]
        }

    def settings = new groovy.json.JsonSlurper().parse(file(params.server_settings))
    def requiredSettings = [
        'db_host', 'db_port', 'db_user', 'db_password', 'db_read_user',
        'registry_host', 'registry_port', 'registry_user',
        'registry_password', 'registry_database'
    ]
    def missingSettings = requiredSettings.findAll { setting -> !settings[setting] }
    if (missingSettings)
        throw new IllegalArgumentException(
            "server_settings is missing required fields: ${missingSettings.join(', ')}"
        )

    db_config_ch = channel
        .value(
            tuple(
                settings.db_host,
                settings.db_port,
                settings.db_user,
                settings.db_password,
                settings.db_read_user
            )
        )

    registry_config_ch = channel
        .value(
            tuple(
                settings.registry_host,
                settings.registry_port,
                settings.registry_user,
                settings.registry_password,
                settings.registry_database
            )
        )

    IMPORT_REFSEQ(
        sample_ch,
        db_config_ch,
        registry_config_ch
    )
}
