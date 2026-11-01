include { VALIDATE_AGGREGATION_MANIFESTS } from '../modules/validate_aggregation_manifests.nf'

workflow PLAN_AGGREGATES {
    take:
    input_manifests // tuple [track_manifest, aggregation_manifest, manifest_base_dir, track_types, allow_opaque_normalization]
    validator

    main:
    planned = VALIDATE_AGGREGATION_MANIFESTS(input_manifests, validator)

    emit:
    plan = planned.plan
    metadata = planned.metadata
    validation = planned.validation
    versions = planned.versions
}
