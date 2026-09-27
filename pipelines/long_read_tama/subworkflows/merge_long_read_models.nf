nextflow.enable.dsl = 2

include { TAMA_MERGE_ACCESSION } from '../modules/tama_merge_accession.nf'
include { TMERGE_MERGE_ACCESSION } from '../modules/tmerge_merge_accession.nf'
include { TAMA_MERGE } from '../modules/tama_merge.nf'
include { TMERGE } from '../modules/tmerge.nf'

workflow MERGE_LONG_READ_MODELS {
    take:
    beds

    main:
    // beds: tuple val(backend), val(meta), val(shard), path(validated_bed)
    accession_inputs = beds.map { backend, meta, _shard, bed -> tuple("${backend}@@${meta.id}", backend, meta.id, bed) }
        .groupTuple()
        .map { _key, backends, accessions, accession_beds -> tuple(backends[0], accessions[0], accession_beds.sort { left, right -> left.name <=> right.name }) }

    if (params.merge_tool == 'tmerge') {
        TMERGE_MERGE_ACCESSION(accession_inputs)
        accession_beds = TMERGE_MERGE_ACCESSION.out.bed
        accession_reports = TMERGE_MERGE_ACCESSION.out.report
        accession_versions = TMERGE_MERGE_ACCESSION.out.versions
    } else {
        TAMA_MERGE_ACCESSION(accession_inputs)
        accession_beds = TAMA_MERGE_ACCESSION.out.bed
        accession_reports = TAMA_MERGE_ACCESSION.out.report
        accession_versions = TAMA_MERGE_ACCESSION.out.versions
    }

    cohort_inputs = accession_beds.map { backend, _accession, bed -> tuple("${backend}@@${params.cohort_id}", backend, params.cohort_id, bed) }
        .groupTuple()
        .map { _key, backends, cohorts, accession_models -> tuple(backends[0], cohorts[0], accession_models.sort { left, right -> left.name <=> right.name }) }
    if (params.merge_tool == 'tmerge') {
        TMERGE(cohort_inputs)
        merged = TMERGE.out.bed
        reports = accession_reports.mix(TMERGE.out.merge_report).mix(TMERGE.out.trans_report)
        versions = accession_versions.mix(TMERGE.out.versions)
    } else {
        TAMA_MERGE(cohort_inputs)
        merged = TAMA_MERGE.out.bed
        reports = accession_reports.mix(TAMA_MERGE.out.gene_report).mix(TAMA_MERGE.out.merge_report).mix(TAMA_MERGE.out.trans_report)
        versions = accession_versions.mix(TAMA_MERGE.out.versions)
    }

    emit:
    bed = merged
    accession_report = accession_reports
    reports
    versions
}
