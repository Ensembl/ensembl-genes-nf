include { VALIDATE_SUPPORTING_EVIDENCE } from '../modules/validate_supporting_evidence.nf'

workflow PLAN_SUPPORTING_EVIDENCE {
    take:
    input_manifest
    validator

    main:
    planned = VALIDATE_SUPPORTING_EVIDENCE(input_manifest, validator)

    emit:
    plan = planned.plan
    metadata = planned.metadata
    validation = planned.validation
    versions = planned.versions
}
