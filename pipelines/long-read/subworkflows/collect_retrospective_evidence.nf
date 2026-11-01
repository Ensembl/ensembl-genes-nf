include { COLLECT_RETROSPECTIVE_EVIDENCE as COLLECT_RETROSPECTIVE_EVIDENCE_PROCESS } from '../modules/retrospective_evidence.nf'

workflow COLLECT_RETROSPECTIVE_EVIDENCE {
    take:
    normalized_manifest
    reference
    cohort_id

    main:
    COLLECT_RETROSPECTIVE_EVIDENCE_PROCESS(normalized_manifest, reference, cohort_id)

    emit:
    evidence = COLLECT_RETROSPECTIVE_EVIDENCE_PROCESS.out.evidence
    reference_identity = COLLECT_RETROSPECTIVE_EVIDENCE_PROCESS.out.reference_identity
    versions = COLLECT_RETROSPECTIVE_EVIDENCE_PROCESS.out.versions
}
