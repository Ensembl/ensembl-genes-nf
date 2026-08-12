#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.sources_csv = null
params.outdir = './chm13_annotation_release_results'
params.assembly_name = 'assembly'
params.run_external_validation = true

include { ANNOTATION_REVIEW } from './subworkflows/annotation_review.nf'

workflow {
    if (!params.sources_csv) {
        error 'Provide --sources_csv with columns: sample_id,projected_gff,manual_gff,decision_tsv,assembly_fasta,assembly_fai'
    }

    integration_script = file("${baseDir}/bin/annotation_integration.py")
    sanitize_script = file("${baseDir}/bin/sanitize_gff.py")

    def sources_csv_file = file(params.sources_csv, checkIfExists: true)
    def sources_base = sources_csv_file.parent

    channel
        .fromPath(sources_csv_file)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample_id]
            def emptyProjected = file("${baseDir}/assets/empty_projected.gff3")
            def emptyManual = file("${baseDir}/assets/empty_manual.gff3")
            def emptyDecisions = file("${baseDir}/assets/empty_decisions.tsv")
            tuple(meta,
                row.projected_gff && row.projected_gff != '-' ? file(row.projected_gff.startsWith('/') ? row.projected_gff : "${sources_base}/${row.projected_gff}", checkIfExists: true) : emptyProjected,
                row.manual_gff && row.manual_gff != '-' ? file(row.manual_gff.startsWith('/') ? row.manual_gff : "${sources_base}/${row.manual_gff}", checkIfExists: true) : emptyManual,
                row.decision_tsv && row.decision_tsv != '-' ? file(row.decision_tsv.startsWith('/') ? row.decision_tsv : "${sources_base}/${row.decision_tsv}", checkIfExists: true) : emptyDecisions,
                file(row.assembly_fasta.startsWith('/') ? row.assembly_fasta : "${sources_base}/${row.assembly_fasta}", checkIfExists: true),
                file(row.assembly_fai.startsWith('/') ? row.assembly_fai : "${sources_base}/${row.assembly_fai}", checkIfExists: true))
        }
        .set { sources }

    ANNOTATION_REVIEW(sources, integration_script, sanitize_script)
}
