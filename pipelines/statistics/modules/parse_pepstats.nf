process PARSE_PEPSTATS {
    label 'process_light'

    tag "${meta.dbname}:pepstats-parse"

    input:
    tuple val(meta), path(pepstats_file)

    output:
    path "versions.yml", emit: versions

    script:
    def pepstats_script = "${projectDir}/bin/pepstats_parser.py"

    """
    python ${pepstats_script} \\
        --pepstats_file ${pepstats_file} \\
        --db_host ${params.host} \\
        --db_port ${params.port} \\
        --db_user ${params.user} \\
        --db_password ${params.password} \\
        --db_name ${meta.dbname}


    PYTHON_VERSION=\$(python --version 2>&1 | awk '{print \$2}')
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$PYTHON_VERSION
    END_VERSIONS
    """

    stub:
    """
    touch versions.yml
    """
}
