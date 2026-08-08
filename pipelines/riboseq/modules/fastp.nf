/*
 * FASTP - Adapter trimming with multiple modes
 *
 * Supports three adapter detection modes:
 * 1. Auto-detection: fastp's built-in adapter detection (no external input needed)
 * 2. Explicit sequence: User-provided adapter sequence via params
 * 3. Adapter FASTA: User-provided or pipeline-generated adapter file
 *
 * The module determines the mode based on inputs and params.
 */

process FASTP {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::fastp=0.23.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--h5f740d0_0' :
        'biocontainers/fastp:0.23.4--h5f740d0_0' }"

    input:
    tuple val(meta), path(reads)
    path adapter_fasta  // Optional: can be empty/placeholder file

    output:
    tuple val(meta), path("*_trimmed.fastq.gz"), emit: trimmed_fastq
    tuple val(meta), path("*_fastp.json"), emit: json
    tuple val(meta), path("*_fastp.html"), emit: html
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_len = params.min_read_length ?: 20
    def max_len_arg = params.max_read_length && params.max_read_length > 0 ? "--length_limit ${params.max_read_length}" : ''

    // Determine adapter arguments based on available inputs and params
    // Priority: 1) explicit sequence param, 2) adapter FASTA file, 3) auto-detection
    def adapter_args = ''
    if (params.fastp_adapter_sequence) {
        // Use explicit adapter sequence from params
        adapter_args = "--adapter_sequence ${params.fastp_adapter_sequence}"
    } else if (adapter_fasta && adapter_fasta.name != 'NO_ADAPTER_FILE' && adapter_fasta.size() > 0) {
        // Use provided adapter FASTA file
        adapter_args = "--adapter_fasta ${adapter_fasta}"
    }
    // If neither, fastp uses auto-detection by default

    """
    fastp \\
        -i ${reads} \\
        -o ${prefix}_trimmed.fastq.gz \\
        ${adapter_args} \\
        --length_required ${min_len} \\
        ${max_len_arg} \\
        --json ${prefix}_fastp.json \\
        --html ${prefix}_fastp.html \\
        --thread ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_trimmed.fastq.gz
    touch ${prefix}_fastp.json
    touch ${prefix}_fastp.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: 0.23.4
    END_VERSIONS
    """
}
