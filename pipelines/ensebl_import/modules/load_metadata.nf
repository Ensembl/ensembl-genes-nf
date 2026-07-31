process LOAD_METADATA {

    tag "${meta.id}"

    input:
    tuple val(meta),
          val(db_name),
          path(sql_file)

    tuple val(db_host),
          val(db_port),
          val(db_user),
          val(db_password)

    script:
    """
    mysql \
        --host ${db_host} \
        --port ${db_port} \
        --user ${db_user} \
        --password=${db_password} \
        ${db_name} < ${sql_file}
    """
}