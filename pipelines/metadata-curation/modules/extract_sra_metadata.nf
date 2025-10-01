process EXTRACT_SRA_METADATA {
    tag "$query"
    label 'process_medium'

    // No conda/container - use base environment

    input:
    val query
    val max_results
    val email

    output:
    path "sra_metadata.json", emit: metadata
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def email_arg = email ? "--email ${email}" : ''
    """
    python ${projectDir}/bin/extract_sra_metadata.py \\
        "${query}" \\
        -o sra_metadata.json \\
        --max-results ${max_results} \\
        ${email_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        requests: \$(python -c "import requests; print(requests.__version__)")
    END_VERSIONS
    """

    stub:
    """
    echo '{"query": "${query}", "bioproject_count": 0, "studies": []}' > sra_metadata.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        requests: \$(python -c "import requests; print(requests.__version__)")
    END_VERSIONS
    """
}