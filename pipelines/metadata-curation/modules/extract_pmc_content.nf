process EXTRACT_PMC_CONTENT {
    tag "$query"
    label 'process_low'

    // No conda/container - use base environment

    input:
    val query
    val max_results
    val email
    val include_full_text

    output:
    path "pmc_content.json", emit: content
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def email_arg = email ? "--email ${email}" : ''
    def full_text_arg = include_full_text ? "--full-text" : ''
    """
    python ${projectDir}/bin/extract_pmc_content.py \\
        "${query}" \\
        -o pmc_content.json \\
        --max-results ${max_results} \\
        ${email_arg} \\
        ${full_text_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        requests: \$(python -c "import requests; print(requests.__version__)")
    END_VERSIONS
    """

    stub:
    """
    echo '{"query": "${query}", "article_count": 0, "articles": []}' > pmc_content.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        requests: \$(python -c "import requests; print(requests.__version__)")
    END_VERSIONS
    """
}