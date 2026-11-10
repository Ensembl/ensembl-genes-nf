process GET_CORE_METADATA {
    label 'python'

    tag "${meta.id}"

    input:
    val meta

    tuple val(db_host),
          val(db_port),
          val(db_user),
          val(db_password),
          val(db_read_user)

    output:
    tuple val(meta), val(db_host), path('core_metadata.json'), emit: metadata
    path 'versions.yml', emit: versions

    script:
    """
    get_core_metadata.py \\
        --database ${meta.db_name} \\
        --host ${db_host} \\
        --port ${db_port} \\
        --user ${db_read_user} \\
        --password '' \\
        > core_metadata.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print \$2}')
        get_core_metadata.py: unknown
    END_VERSIONS
    """

    stub:
    """
    printf '%s\n' '{"genebuild_version":"EXT01","last_geneset_update":"2026-01"}' > core_metadata.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
        get_core_metadata.py: stub
    END_VERSIONS
    """
}
