process DETECT_BIOGROUPS {
    tag "$meta.id"
    label 'process_low'

    // No conda/container - use base environment

    input:
    tuple val(meta), path(processed_metadata)

    output:
    tuple val(meta), path("biogroups.json"), emit: biogroups
    path "biogroup_report.json"            , emit: report
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    python ${projectDir}/bin/biogroup_detector.py \\
        ${processed_metadata} \\
        -o biogroups.json \\
        --report biogroup_report.json \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    echo '{"biogroup_count": 0, "total_samples": 0, "biogroups": []}' > biogroups.json
    echo '{"total_biogroups": 0}' > biogroup_report.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}