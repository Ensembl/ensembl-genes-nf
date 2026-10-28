process LOAD_TAXONOMY {
    label 'python'

    tag "${meta.id}"

    publishDir { "${params.outdir}/metadata/${meta.id}" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(loaded_marker)
    tuple val(db_host),
          val(db_port),
          val(db_user),
          val(db_password),
          val(db_read_user)

    output:
    tuple val(meta), path("${meta.id}.taxonomy.done"), emit: taxonomy
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    get_taxonomy.py \\
        --database ${meta.db_name} \\
        --host ${db_host} \\
        --port '${db_port}' \\
        --user_w ${db_user} \\
        --password '${db_password}' \\
        --read_user ${db_read_user} \\
        --core_mode

    touch ${meta.id}.taxonomy.done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print \$2}')
        get_taxonomy.py: unknown
    END_VERSIONS

    """

    stub:
    """
    touch ${meta.id}.taxonomy.done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
        get_taxonomy.py: stub
    END_VERSIONS
    """
}
