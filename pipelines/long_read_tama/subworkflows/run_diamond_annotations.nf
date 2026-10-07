nextflow.enable.dsl = 2

include { GTF_TO_BED12 } from '../modules/gtf_to_bed12.nf'
include { RUN_DIAMOND_QC } from './run_diamond_qc.nf'

workflow RUN_DIAMOND_ANNOTATIONS {
    take:
    annotation_manifest
    reference

    main:
    rows = annotation_manifest
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def annotation = row.annotation ?: row.path
            def backend = row.backend?.toString()?.trim()
            def scope = row.scope?.toString()?.trim()
            def identifier = row.id ?: row.sample
            def format = (row.format ?: 'gtf').toString().toLowerCase()
            if (!annotation || !identifier || !backend || !scope)
                error 'diamond_annotation_manifest requires id, backend, scope, annotation, and optional format columns'
            if (!(format in ['gtf', 'gff', 'gff3', 'bed12']))
                error "Unsupported annotation format '${format}' for ${identifier}/${backend}/${scope}"
            def meta = [id: identifier.toString(), backend: backend, scope: scope]
            tuple(meta, backend, scope, format, file(annotation, checkIfExists: true))
        }

    gtf_inputs = rows
        .filter { _meta, _backend, _scope, format, _annotation -> format in ['gtf', 'gff', 'gff3'] }
        .map { meta, backend, scope, _format, annotation -> tuple(meta, backend, scope, annotation) }
    bed_inputs = rows
        .filter { _meta, _backend, _scope, format, _annotation -> format == 'bed12' }
        .map { meta, backend, _scope, _format, annotation -> tuple(meta, backend, annotation) }

    GTF_TO_BED12(gtf_inputs)
    converted_beds = GTF_TO_BED12.out.bed
        .map { meta, backend, _scope, bed -> tuple(meta, backend, bed) }
    annotation_beds = converted_beds.mix(bed_inputs)

    RUN_DIAMOND_QC(annotation_beds, reference)

    emit:
    report = RUN_DIAMOND_QC.out.report
    summary = RUN_DIAMOND_QC.out.summary
    provenance = RUN_DIAMOND_QC.out.provenance
    versions = GTF_TO_BED12.out.versions.mix(RUN_DIAMOND_QC.out.versions)
}
