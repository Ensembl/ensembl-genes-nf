process VALIDATE_AGGREGATION_MANIFESTS {
    tag 'aggregation-plan'
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    tuple path(track_manifest), path(aggregation_manifest), val(manifest_base_dir), val(track_types), val(allow_opaque_normalization)
    path validator

    output:
    path 'aggregation_plan.tsv', emit: plan
    path 'aggregate_metadata.json', emit: metadata
    path 'aggregation_validation.tsv', emit: validation
    path 'aggregation_plan.versions.yml', emit: versions

    script:
    def opaque_flag = allow_opaque_normalization ? '--allow-opaque-normalization' : ''
    """
    python3 ${validator} \
        --track-manifest ${track_manifest} \
        --aggregation-manifest ${aggregation_manifest} \
        --manifest-base-dir '${manifest_base_dir}' \
        --track-types '${track_types}' \
        ${opaque_flag} \
        --output-plan aggregation_plan.tsv \
        --output-metadata aggregate_metadata.json \
        --output-validation aggregation_validation.tsv
    printf '"%s":\\n    validate_aggregation_manifests: pipeline\\n' '${task.process}' > aggregation_plan.versions.yml
    """

    stub:
    def opaque_flag = allow_opaque_normalization ? '--allow-opaque-normalization' : ''
    """
    python3 ${validator} \\
        --track-manifest ${track_manifest} \\
        --aggregation-manifest ${aggregation_manifest} \\
        --manifest-base-dir '${manifest_base_dir}' \\
        --track-types '${track_types}' \\
        ${opaque_flag} \\
        --output-plan aggregation_plan.tsv \\
        --output-metadata aggregate_metadata.json \\
        --output-validation aggregation_validation.tsv
    printf '"%s":\\n    validate_aggregation_manifests: stub\\n' '${task.process}' > aggregation_plan.versions.yml
    """
}
