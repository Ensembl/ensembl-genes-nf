#!/usr/bin/env nextflow

include { validateParameters; samplesheetToList } from 'plugin/nf-schema'
include { QC } from './workflows/qc.nf'
include { PAIRWISE_ANNOTATION_COMPARISON } from './subworkflows/pairwise_annotation/pairwise_annotation_comparison.nf'

workflow {
    // nf-schema validates the parameter and annotation manifest contracts.
    validateParameters()

    if (params.mode == 'pairwise') {
        pair_csv = params.pairwise_csv ?: params.input_csv
        channel
            .fromPath(pair_csv, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def pairId = row.pair_id ?: row.sample ?: row.id
                def assembly = row.assembly ?: row.assembly_accession ?: pairId
                def sourceA = row.source_a ?: 'source_a'
                def sourceB = row.source_b ?: 'source_b'
                def sourceAPath = row.source_a_path ?: row.gff_a ?: row.annotation_a
                def sourceBPath = row.source_b_path ?: row.gff_b ?: row.annotation_b
                def formatA = row.source_a_format ?: row.format_a ?: 'auto'
                def formatB = row.source_b_format ?: row.format_b ?: 'auto'
                if (formatA.toLowerCase() == 'auto')
                    formatA = sourceAPath.toLowerCase().contains('.bed') ? 'bed12' : (sourceAPath.toLowerCase().contains('.gtf') ? 'gtf' : 'gff3')
                if (formatB.toLowerCase() == 'auto')
                    formatB = sourceBPath.toLowerCase().contains('.bed') ? 'bed12' : (sourceBPath.toLowerCase().contains('.gtf') ? 'gtf' : 'gff3')
                def report = row.assembly_report ?: ''
                def evidenceA = row.evidence_a ?: ''
                def evidenceB = row.evidence_b ?: ''
                def referenceAnnotation = row.reference_annotation ?: row.reference_path ?: params.pairwise_reference_annotation ?: ''
                def referenceFormat = (row.reference_format ?: params.pairwise_reference_format ?: 'auto').toLowerCase()
                [sourceAPath, sourceBPath, report, evidenceA, evidenceB, referenceAnnotation].eachWithIndex { value, index ->
                    if (value && !value.toString().matches('(?i)^(https?|ftp)://.*') && !value.toString().startsWith('/')) {
                        def resolved = "${projectDir}/${value}"
                        if (index == 0) sourceAPath = resolved
                        if (index == 1) sourceBPath = resolved
                        if (index == 2) report = resolved
                        if (index == 3) evidenceA = resolved
                        if (index == 4) evidenceB = resolved
                        if (index == 5) referenceAnnotation = resolved
                    }
                }
                if (!pairId || !sourceAPath || !sourceBPath)
                    error "Pairwise CSV requires pair_id, source_a_path and source_b_path (or legacy aliases)"
                def meta = [
                    id: pairId, pair_id: pairId, assembly: assembly,
                    source_a: sourceA, source_b: sourceB,
                    format_a: formatA.toLowerCase(), format_b: formatB.toLowerCase(),
                    source_a_assembly: row.source_a_assembly ?: '', source_b_assembly: row.source_b_assembly ?: '',
                    evidence_a: evidenceA, evidence_b: evidenceB,
                    reference_annotation: referenceAnnotation, reference_format: referenceFormat
                ]
                tuple(meta, sourceAPath, sourceBPath, report, evidenceA, evidenceB, referenceAnnotation)
            }
            .set { pairwise_samples_ch }
        PAIRWISE_ANNOTATION_COMPARISON(pairwise_samples_ch)
    }
    else {

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
}
