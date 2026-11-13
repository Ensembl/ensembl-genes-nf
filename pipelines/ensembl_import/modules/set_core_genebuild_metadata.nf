process SET_CORE_GENEBUILD_METADATA {
    label 'process_light'

    tag "${meta.id}"

    input:
    tuple val(meta), path(loaded_marker)

    tuple val(db_host),
          val(db_port),
          val(db_user),
          val(db_password)

    output:
    tuple val(meta), val(db_host), path('core_genebuild_metadata.json'), emit: updated
    path 'versions.yml', emit: versions

    script:
    """
    set_core_genebuild_metadata.py \\
        --database '${meta.db_name}' \\
        --host '${db_host}' \\
        --port '${db_port}' \\
        --user '${db_user}' \\
        --password '${db_password}' \\
        --genebuild-version EXT01 \\
        > core_genebuild_metadata.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print \$2}')
        set_core_genebuild_metadata.py: unknown
    END_VERSIONS
    """

    stub:
    """
    printf '%s\n' '{"genebuild_version":"EXT01","last_geneset_update":"2026-01"}' > core_genebuild_metadata.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
        set_core_genebuild_metadata.py: stub
    END_VERSIONS
    """
}
