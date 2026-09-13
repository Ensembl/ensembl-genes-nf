process GET_SAMPLE_GENE {
    label 'process_light'

    tag "${meta.id}"

    input:
    tuple val(meta), path(loaded_marker)

    tuple val(db_host),
          val(db_port),
          val(db_user),
          val(db_password)

    output:
    tuple val(meta), path("${meta.id}.sample_gene.done"), emit: sample_gene
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def script = "${projectDir}/bin/get_sample_gene.py"

    """
    python ${script} \\
        --db-name ${meta.db_name} \\
        --host ${db_host} \\
        --port ${db_port} \\
        --user ${db_user} \\
        --password '${db_password}' \\
        ${args}

    touch ${meta.id}.sample_gene.done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print \$2}')
        get_sample_gene.py: unknown
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sample_gene.done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
        get_sample_gene.py: stub
    END_VERSIONS
    """
}
