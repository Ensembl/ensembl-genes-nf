include { VALIDATE_INPUT_MANIFEST } from '../modules/preparation.nf'

workflow VALIDATE_RETROSPECTIVE_INPUT {
    take:
    manifest
    mode
    cohort_id
    manifest_root
    ccs_mode

    main:
    VALIDATE_INPUT_MANIFEST(manifest, mode, cohort_id, manifest_root, ccs_mode, false, '')

    emit:
    normalized = VALIDATE_INPUT_MANIFEST.out.normalized
    versions = channel.empty()
}
