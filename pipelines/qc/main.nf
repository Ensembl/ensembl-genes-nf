#!/usr/bin/env nextflow

include { validateParameters; samplesheetToList } from 'plugin/nf-schema'
include { QC } from './workflows/qc.nf'

workflow {
    // nf-schema validates the parameter and annotation manifest contracts.
    validateParameters()

    // This is pipeline-owned configuration, not a user-supplied input.
    feature_levels_yaml = file("${projectDir}/assets/feature_levels.yaml")

    def samplesheet_schema = [
        agat:          'assets/agat_samplesheet.json',
        interproscan:  'assets/protein_qc_samplesheet.json',
        diamond:       'assets/full_annotation_samplesheet.json',
        combined:      'assets/full_annotation_samplesheet.json'
    ][params.mode]

    annotation_ch = channel.fromList(
        samplesheetToList(params.input_csv, samplesheet_schema)
    )

    QC(
        annotation_ch,
        params.mode,
        feature_levels_yaml,
        params.database,
        params.data_file_path,
        params.diamond_reference_proteins,
        params.diamond_reference_db
    )
}
