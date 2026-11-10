include { PREPARE_RETROSPECTIVE_EVIDENCE; PACKAGE_RETROSPECTIVE_EVIDENCE; VERIFY_RETROSPECTIVE_EVIDENCE } from '../modules/retrospective_evidence.nf'

workflow COLLECT_RETROSPECTIVE_EVIDENCE {
    take:
    normalized_manifest
    reference
    cohort_id

    main:
    PREPARE_RETROSPECTIVE_EVIDENCE(normalized_manifest, reference, cohort_id)
    PACKAGE_RETROSPECTIVE_EVIDENCE(
        PREPARE_RETROSPECTIVE_EVIDENCE.out.models,
        PREPARE_RETROSPECTIVE_EVIDENCE.out.alignment_manifests,
        normalized_manifest,
        PREPARE_RETROSPECTIVE_EVIDENCE.out.reference_identity,
        PREPARE_RETROSPECTIVE_EVIDENCE.out.model_records,
        cohort_id
    )
    VERIFY_RETROSPECTIVE_EVIDENCE(PACKAGE_RETROSPECTIVE_EVIDENCE.out.evidence, cohort_id)

    emit:
    evidence = VERIFY_RETROSPECTIVE_EVIDENCE.out.evidence
    reference_identity = PREPARE_RETROSPECTIVE_EVIDENCE.out.reference_identity
    versions = channel.empty()
}
