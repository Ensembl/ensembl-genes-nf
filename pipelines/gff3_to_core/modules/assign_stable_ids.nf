process ASSIGN_STABLE_IDS {
    label 'process_single'

    conda "conda-forge::python=3.11 conda-forge::pymysql=1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    input:
    path stats  // load_stats.json from LOAD_GFF3_TO_CORE (used as ordering dependency)

    output:
    path 'stable_id_stats.json', emit: stats
    path 'versions.yml',         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    assign_stable_ids.py \\
        --host       ${params.db_host} \\
        --port       ${params.db_port} \\
        --user       ${params.db_user} \\
        --password   ${params.db_password} \\
        --dbname     ${params.db_name} \\
        --prefix     '${params.stable_id_prefix}' \\
        --species-id ${params.species_id}

    echo '{"stable_ids_assigned":true}' > stable_id_stats.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '{"stable_ids_assigned":true}' > stable_id_stats.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
