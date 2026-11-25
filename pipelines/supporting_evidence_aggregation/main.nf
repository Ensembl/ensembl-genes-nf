#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.enable.strict = true

include { validateParameters } from 'plugin/nf-schema'
include { PLAN_SUPPORTING_EVIDENCE } from './subworkflows/plan_inputs.nf'
include { BUILD_COVERAGE_MATRICES } from './subworkflows/coverage.nf'
include { BUILD_SPLICE_MATRIX_WORKFLOW } from './subworkflows/splice.nf'
include { WRITE_SUPPORTING_EVIDENCE_REPORTS_WORKFLOW } from './subworkflows/reporting.nf'

workflow {
    validateParameters()
    if (!params.track_manifest) error '--track_manifest is required'

    def track_manifest = file(params.track_manifest, checkIfExists: true)
    def manifest_base_dir = params.manifest_base_dir ?
        file(params.manifest_base_dir, checkIfExists: true).toString() :
        track_manifest.toFile().absoluteFile.parentFile.toString()

    input_manifest = channel.of(tuple(
        track_manifest,
        manifest_base_dir,
        params.track_types,
        params.source_entity_type,
        params.species_name ?: '',
        params.allow_opaque_normalization,
    ))
    plan = PLAN_SUPPORTING_EVIDENCE(
        input_manifest,
        file("${projectDir}/bin/validate_supporting_evidence_inputs.py", checkIfExists: true),
    )

    coverage = BUILD_COVERAGE_MATRICES(plan.plan)
    splice = BUILD_SPLICE_MATRIX_WORKFLOW(
        plan.plan,
        file("${projectDir}/bin/build_splice_matrix.py", checkIfExists: true),
    )

    result_files = coverage.outputs.mix(splice.outputs).collect()
    WRITE_SUPPORTING_EVIDENCE_REPORTS_WORKFLOW(
        plan.plan,
        result_files,
        file("${projectDir}/bin/write_supporting_evidence_reports.py", checkIfExists: true),
    )
}
