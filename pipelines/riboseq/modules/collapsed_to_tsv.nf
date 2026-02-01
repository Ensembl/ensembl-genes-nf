/*
 * Convert collapsed FASTA to sorted TSV
 * Simple transformation step for downstream matrix building
 */

process COLLAPSED_TO_TSV {
    tag "${meta.id}"
    label 'process_low'

    // Pure Python - no external dependencies
    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(collapsed_fasta)

    output:
    tuple val(meta), path("${meta.id}.tsv"), emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    collapsed_to_tsv.py \\
        ${collapsed_fasta} \\
        --output ${meta.id}.tsv \\
        --sample-id ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: unknown
    END_VERSIONS
    """
}
