process MAP_ONTOLOGY_TERMS {
    tag "$meta.id"
    label 'process_medium'

    // No conda/container - use base environment

    input:
    tuple val(meta), path(metadata_json)
    val ontology_cache_dir

    output:
    tuple val(meta), path("mapped_metadata.json"), emit: mapped_metadata
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def cache_arg = ontology_cache_dir ? "--cache-dir ${ontology_cache_dir}" : ''
    """
    python ${projectDir}/bin/ontology_mapper.py \\
        ${metadata_json} \\
        -o mapped_metadata.json \\
        ${cache_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        requests: \$(python -c "import requests; print(requests.__version__)")
    END_VERSIONS
    """

    stub:
    """
    cp ${metadata_json} mapped_metadata.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        requests: \$(python -c "import requests; print(requests.__version__)")
    END_VERSIONS
    """
}