include { FINALIZE_RUN_MANIFEST; VERIFY_PUBLICATION_CHECKSUMS; COLLECT_SOFTWARE_VERSIONS } from '../modules/finalization.nf'

workflow FINALIZE_LONG_READ_PIPELINE {
    take:
    normalized_pair
    shard_plans
    model_reports
    reference_fasta
    reference_identity
    versions
    run_configuration
    stage_status
    evidence_status
    model_publications
    merge_provenance

    main:
    FINALIZE_RUN_MANIFEST(normalized_pair, shard_plans, model_reports, reference_fasta, reference_identity, run_configuration, stage_status, evidence_status, model_publications, merge_provenance)
    VERIFY_PUBLICATION_CHECKSUMS(normalized_pair, shard_plans, model_reports, reference_fasta, reference_identity, FINALIZE_RUN_MANIFEST.out.checksums)
    versions_text = versions.map { version_file -> version_file.text }.collect().map { records -> records.sort().join('\n') }
    COLLECT_SOFTWARE_VERSIONS(versions_text)
    emit:
    run_manifest = FINALIZE_RUN_MANIFEST.out.run_manifest
    checksums = FINALIZE_RUN_MANIFEST.out.checksums
    checksum_verification = VERIFY_PUBLICATION_CHECKSUMS.out.checksum_verification
    software_versions = COLLECT_SOFTWARE_VERSIONS.out.software_versions
}
