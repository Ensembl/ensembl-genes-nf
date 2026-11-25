process VALIDATE_SUPPORTING_EVIDENCE {
    tag 'supporting-evidence-plan'
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'
    publishDir "${params.outdir}/plan", mode: 'copy', overwrite: true

    input:
    tuple path(track_manifest), val(manifest_base_dir), val(track_types), val(entity_type), val(species_name), val(allow_opaque_normalization)
    path validator

    output:
    path 'supporting_evidence_plan.tsv', emit: plan
    path 'supporting_evidence_metadata.json', emit: metadata
    path 'supporting_evidence_validation.tsv', emit: validation
    path 'supporting_evidence_plan.versions.yml', emit: versions

    script:
    def opaque_flag = allow_opaque_normalization ? '--allow-opaque-normalization' : ''
    """
    python3 ${validator} \
        --track-manifest ${track_manifest} \
        --manifest-base-dir '${manifest_base_dir}' \
        --track-types '${track_types}' \
        --entity-type '${entity_type}' \
        --species-name '${species_name}' \
        ${opaque_flag} \
        --output-plan supporting_evidence_plan.tsv \
        --output-metadata supporting_evidence_metadata.json \
        --output-validation supporting_evidence_validation.tsv
    printf '"%s":\n    validate_supporting_evidence: python\n' '${task.process}' > supporting_evidence_plan.versions.yml
    """

    stub:
    def opaque_flag = allow_opaque_normalization ? '--allow-opaque-normalization' : ''
    """
    python3 ${validator} \
        --track-manifest ${track_manifest} \
        --manifest-base-dir '${manifest_base_dir}' \
        --track-types '${track_types}' \
        --entity-type '${entity_type}' \
        --species-name '${species_name}' \
        ${opaque_flag} \
        --output-plan supporting_evidence_plan.tsv \
        --output-metadata supporting_evidence_metadata.json \
        --output-validation supporting_evidence_validation.tsv
    printf '"%s":\n    validate_supporting_evidence: stub\n' '${task.process}' > supporting_evidence_plan.versions.yml
    """
}
