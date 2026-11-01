include { WRITE_MODEL_RECORDS; PACKAGE_LONG_READ_EVIDENCE; VERIFY_LONG_READ_EVIDENCE } from '../modules/evidence.nf'

workflow COLLECT_LONG_READ_EVIDENCE_PRODUCTS {
    take:
    validated_models
    normalized_manifest
    alignment_manifests
    reference_identity
    cohort_id

    main:
    model_records = validated_models.map { meta, backend, scope, bed12, model_manifest, validation ->
        [bed12_name: bed12.name, backend: backend, scope: scope, run_id: meta.run_id, cohort_id: meta.cohort_id, model_manifest: model_manifest.text, validation_report: validation.text]
    }.collect()
    model_records_json = model_records.map { groovy.json.JsonOutput.toJson(it) }
    beds = validated_models.map { _meta, _backend, _scope, bed12, _manifest, _validation -> bed12 }.collect()
    manifest = normalized_manifest.map { normalized, _checksum, _input_checksums, _status -> normalized }
    alignment_files = alignment_manifests
    reference = reference_identity
    WRITE_MODEL_RECORDS(model_records_json, cohort_id)
    PACKAGE_LONG_READ_EVIDENCE(beds, manifest, alignment_files, reference, WRITE_MODEL_RECORDS.out.records, cohort_id)
    VERIFY_LONG_READ_EVIDENCE(PACKAGE_LONG_READ_EVIDENCE.out.evidence, cohort_id)

    emit:
    evidence = VERIFY_LONG_READ_EVIDENCE.out.evidence
    versions = WRITE_MODEL_RECORDS.out.versions.mix(PACKAGE_LONG_READ_EVIDENCE.out.versions).mix(VERIFY_LONG_READ_EVIDENCE.out.versions)
}
