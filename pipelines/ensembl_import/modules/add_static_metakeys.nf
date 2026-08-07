process ADD_STATIC_METAKEYS {
    label 'process_light'

    tag "${meta.id}"
    publishDir "${params.outdir}/metadata/${meta.id}", mode: 'copy', overwrite: true

    input:
    val(meta)

    tuple val(db_host),
          val(db_port),
          val(db_user),
          val(db_password)

    output:
    tuple val(meta), path("${meta.id}.add_static_metakeys.done"), emit: loaded
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def script = "${projectDir}/bin/add_static_metakeys.py"
    def static_json = "${projectDir}/config/static_metakeys.json"

    """
    python ${script} \\
        --db_name ${meta.db_name} \\
        --db_host ${db_host} \\
        --db_port ${db_port} \\
        --db_user ${db_user} \\
        --db_password '${db_password}' \\
        --json_file ${static_json} \\
        --source "refseq"

    touch ${meta.id}.add_static_metakeys.done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print \$2}')
        add_static_metakeys.py: unknown
    END_VERSIONS

    """

    stub:
    """
    touch ${meta.id}.add_static_metakeys.done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
        add_static_metakeys.py: stub
    END_VERSIONS
    """
}
