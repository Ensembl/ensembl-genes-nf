#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.enable.strict = true

include { validateParameters } from 'plugin/nf-schema'
include { COLLECT_RETROSPECTIVE_EVIDENCE } from './subworkflows/collect_retrospective_evidence.nf'
include { VALIDATE_RETROSPECTIVE_INPUT } from './subworkflows/validate_retrospective_input.nf'
include { FINALIZE_LONG_READ_PIPELINE } from './subworkflows/evidence_finalize.nf'

def resolve_reference_path() {
    if (params.reference_fasta) return params.reference_fasta
    if (!params.reference_metadata) error 'retrospective mode requires --reference_fasta or --reference_metadata'
    def metadata_file = file(params.reference_metadata, checkIfExists: true)
    def metadata = new groovy.json.JsonSlurper().parse(metadata_file)
    def path_value = metadata.reference_fasta ?: metadata.fasta ?: metadata.path
    if (!path_value) error '--reference_metadata must contain reference_fasta, fasta, or path'
    path_value.toString().startsWith('/') ? path_value.toString() : file("${metadata_file.parent}/${path_value}").toString()
}

workflow {
    validateParameters(parameters_schema: 'retrospective_schema.json')
    if (!params.retrospective_manifest) error '--retrospective_manifest is required'
    if (params.evidence_mode != 'enabled') error 'retrospective mode requires --evidence_mode enabled'
    if (!params.reference_fasta && !params.reference_metadata) error 'retrospective mode requires --reference_fasta or --reference_metadata'
    manifest = file(params.retrospective_manifest, checkIfExists: true)
    reference = file(resolve_reference_path(), checkIfExists: true)
    if (params.reference_fai) file(params.reference_fai, checkIfExists: true)
    VALIDATE_RETROSPECTIVE_INPUT(manifest, 'retrospective', params.cohort_id, manifest.parent.toString(), params.ccs_mode)
    normalized_manifest = VALIDATE_RETROSPECTIVE_INPUT.out.normalized.map { normalized, _checksum, _input_checksums, _status -> normalized }
    COLLECT_RETROSPECTIVE_EVIDENCE(normalized_manifest, reference, params.cohort_id)
    retrospective_artifact = VALIDATE_RETROSPECTIVE_INPUT.out.normalized.map { _normalized, _checksum, _input_checksums, status -> status }
    configuration = groovy.json.JsonOutput.toJson([pipeline_version: workflow.revision ?: 'uncommitted', nextflow_version: workflow.nextflow.version.toString(), mode: 'retrospective', cohort_id: params.cohort_id, evidence_mode: 'enabled', ccs_mode: params.ccs_mode, retrospective_manifest: manifest.toString(), reference_fasta: reference.toString(), reference_metadata: params.reference_metadata])
    stage_status = groovy.json.JsonOutput.toJson([preparation: 'COMPLETE', alignment: 'NOT_APPLICABLE', sharding: 'NOT_APPLICABLE', reconstruction: 'NOT_APPLICABLE', merge: 'NOT_APPLICABLE', normalization: 'NOT_APPLICABLE', validation: 'COMPLETE', orf: 'DISABLED', diamond: 'DISABLED', evidence: 'COMPLETE', finalization: 'COMPLETE'])
    evidence_status = COLLECT_RETROSPECTIVE_EVIDENCE.out.evidence.map { _evidence -> 'COMPLETE' }
    normalized_pair = VALIDATE_RETROSPECTIVE_INPUT.out.normalized.map { normalized, checksum, input_checksums, status -> tuple(normalized, input_checksums, status) }
    empty_records = channel.of(file("${projectDir}/resources/empty_jsonl.jsonl", checkIfExists: true))
    FINALIZE_LONG_READ_PIPELINE(normalized_pair, retrospective_artifact, retrospective_artifact, reference, COLLECT_RETROSPECTIVE_EVIDENCE.out.reference_identity, VALIDATE_RETROSPECTIVE_INPUT.out.versions.mix(COLLECT_RETROSPECTIVE_EVIDENCE.out.versions), configuration, stage_status, evidence_status, empty_records, empty_records)
}
