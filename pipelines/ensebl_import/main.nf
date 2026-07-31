#!/usr/bin/env nextflow

import groovy.json.JsonSlurper

include { validateParameters } from 'plugin/nf-schema'
include { IMPORT_REFSEQ } from './subworkflows/import_refseq.nf'

def validate_params() {
    def errors = []
    if (!params.input_csv)
        errors << "  --input_csv is required"
    else if (!file(params.input_csv).exists())
        errors << "  --input_csv does not exist: ${params.input_csv}"

    if (!params.server_settings)
        errors << "  --server_settings is required"
    else if (!file(params.server_settings).exists())
        errors << "  --server_settings does not exist: ${params.server_settings}"

    if (!params.ensembl_genes_repo)
        errors << "  --ensembl_genes_repo is required"
    else if (!file(params.ensembl_genes_repo).exists())
        errors << "  --ensembl_genes_repo does not exist: ${params.ensembl_genes_repo}"

    if (errors) {
        log.error "Missing or invalid parameters:\n${errors.join('\n')}"
        System.exit(1)
    }
}

workflow {

    validateParameters()
    validate_params()

    //Read CSV with columns: gcf, species
    Channel
        .fromPath(params.input_csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            // row is a map: [gcf: 'GCF23442512.1', species: 'Homo sapiens']
            def meta = [
                id     : row.gcf,
                species: row.species
            ]
            tuple(meta)
        }
        .set { sample_ch }

    // Read server_settings file

    def settings = new JsonSlurper().parse(file(params.server_settings))

    Channel
        .value(
            tuple(
                settings.db_host,
                settings.db_port,
                settings.db_user,
                settings.db_password
            )
        )
        .set { db_config_ch }

    // Repository location

    ensembl_genes_repo = file(params.ensembl_genes_repo)

    // Run import_refseq pipeline
    IMPORT_REFSEQ(
        sample_ch,
        db_config_ch
    )
}
