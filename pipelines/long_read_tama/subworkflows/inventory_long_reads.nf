nextflow.enable.dsl = 2

include { VALIDATE_LONG_READ_MANIFEST } from '../modules/validate_manifest.nf'
include { DISCOVER_LONG_READ_DATA } from '../modules/discover_long_read.nf'
include { INSPECT_LONG_READ_MANIFEST } from '../modules/inspect_manifest.nf'
include { RESOLVE_LONG_READ_METADATA } from '../modules/resolve_metadata.nf'
include { AUTO_APPROVE_SELECTED_LONG_READS } from '../modules/auto_approve_selected.nf'
include { TRANSCRIPTOMIC_DATA } from '../modules/transcriptomic_data.nf'

workflow INVENTORY_LONG_READS {
    take:
    manifest
    taxon_id
    discovery_script
    classifier_script
    resolver_script
    accession_extractor
    normaliser_script
    approver_script
    validator_script
    metadata_json
    discovery_cache_dir
    metadata_cache_dir
    auto_approve

    main:
    discovered_manifest = channel.empty()
    if (taxon_id) {
        TRANSCRIPTOMIC_DATA(taxon_id)
        DISCOVER_LONG_READ_DATA(
            taxon_id,
            discovery_script,
            classifier_script,
            resolver_script,
            TRANSCRIPTOMIC_DATA.out.candidates,
            discovery_cache_dir
        )
        discovered_manifest = DISCOVER_LONG_READ_DATA.out.results.map { discovery_dir ->
            file("${discovery_dir}/proposed_long_read_manifest.tsv")
        }
    }

    candidate_manifest = taxon_id ? discovered_manifest : manifest
    VALIDATE_LONG_READ_MANIFEST(candidate_manifest, normaliser_script)

    if (metadata_json) {
        inventory_metadata = metadata_json
    } else {
            RESOLVE_LONG_READ_METADATA(
            VALIDATE_LONG_READ_MANIFEST.out.manifest,
            resolver_script,
            metadata_cache_dir,
            accession_extractor
        )
        inventory_metadata = RESOLVE_LONG_READ_METADATA.out.metadata
    }

    INSPECT_LONG_READ_MANIFEST(
        VALIDATE_LONG_READ_MANIFEST.out.manifest,
        inventory_metadata,
        validator_script
    )

    approved_manifest = channel.empty()
    approval_audit = channel.empty()
    version_ch = VALIDATE_LONG_READ_MANIFEST.out.versions
        .mix(INSPECT_LONG_READ_MANIFEST.out.versions)

    if (taxon_id)
        version_ch = TRANSCRIPTOMIC_DATA.out.versions.mix(DISCOVER_LONG_READ_DATA.out.versions).mix(version_ch)

    if (auto_approve) {
        AUTO_APPROVE_SELECTED_LONG_READS(
            INSPECT_LONG_READ_MANIFEST.out.reports,
            approver_script
        )
        approved_manifest = AUTO_APPROVE_SELECTED_LONG_READS.out.manifest
        approval_audit = AUTO_APPROVE_SELECTED_LONG_READS.out.audit
        version_ch = version_ch.mix(AUTO_APPROVE_SELECTED_LONG_READS.out.versions)
    }

    if (!metadata_json)
        version_ch = version_ch.mix(RESOLVE_LONG_READ_METADATA.out.versions)

    emit:
    approved = approved_manifest
    reports = INSPECT_LONG_READ_MANIFEST.out.reports
    audit = approval_audit
    versions = version_ch
}
