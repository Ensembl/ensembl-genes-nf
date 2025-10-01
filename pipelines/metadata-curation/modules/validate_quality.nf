process VALIDATE_QUALITY {
    tag "$meta.id"
    label 'process_low'

    // No conda/container - use base environment

    input:
    tuple val(meta), path(processed_metadata)
    path biogroups

    output:
    path "validation_report.json", emit: report
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def biogroups_arg = biogroups.name != 'NO_FILE' ? "--biogroups ${biogroups}" : ''
    """
    python ${projectDir}/bin/quality_validator.py \\
        ${processed_metadata} \\
        -o validation_report.json \\
        ${biogroups_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    echo '{"validation_summary": {"total_issues": 0, "errors": 0, "warnings": 0}}' > validation_report.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}