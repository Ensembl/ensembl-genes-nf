process EXTRACT_STUDY_CONTEXT {
    tag "$meta.id"
    label 'process_low'

    // No conda/container - use base environment

    input:
    tuple val(meta), path(mapped_metadata)

    output:
    tuple val(meta), path("context_metadata.json"), emit: context_metadata
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    python ${projectDir}/bin/llm_context_extractor.py \\
        ${mapped_metadata} \\
        -o context_metadata.json \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    cp ${mapped_metadata} context_metadata.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}