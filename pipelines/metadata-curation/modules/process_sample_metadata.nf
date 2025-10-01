process PROCESS_SAMPLE_METADATA {
    tag "$meta.id"
    label 'process_medium'

    // No conda/container - use base environment

    input:
    tuple val(meta), path(context_metadata)

    output:
    tuple val(meta), path("processed_metadata.json"), emit: processed_metadata
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    python ${projectDir}/bin/evidence_metadata_processor.py \\
        ${context_metadata} \\
        -o processed_metadata.json \\
        --field-config ${projectDir}/conf/field_registry.json \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    cp ${context_metadata} processed_metadata.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}