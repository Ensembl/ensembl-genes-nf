process LOAD_METADATA {
    label 'process_light'

    tag "${meta.id}"

    input:
    tuple val(meta),
          path(sql_file)

    tuple val(db_host),
          val(db_port),
          val(db_user),
          val(db_password)

    output:
    tuple val(meta), path("${meta.id}.load_metadata.done"), emit: loaded
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mysql \
        --host ${db_host} \
        --port ${db_port} \
        --user ${db_user} \
        --password='${db_password}' \
        ${meta.db_name} < ${sql_file}

    touch ${meta.id}.load_metadata.done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mysql: \$(mysql --version | sed 's/^.*Distrib //' | awk '{print \$1}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.load_metadata.done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mysql: stub
    END_VERSIONS
    """
}
