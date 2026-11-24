/*
 * CUTADAPT_RPF - Traditional single-end Ribo-Seq adapter trimming.
 *
 * Cutadapt is the explicit multi-adapter alternative.  A FASTA adapter file
 * is passed through cutadapt's file: syntax, preserving all records.
 */

process CUTADAPT_RPF {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::cutadapt=5.2"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:5.2--py39hbcbf7aa_0' :
        'biocontainers/cutadapt:5.2--py39hbcbf7aa_0' }"

    input:
    tuple val(meta), path(reads), path(adapter_fasta)

    output:
    tuple val(meta), path("*_trimmed.fastq.gz"), emit: trimmed_fastq
    tuple val(meta), path("*_cutadapt.log"), emit: log
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_len = params.min_read_length ?: 20
    def max_len = params.max_read_length ?: 0
    def max_len_arg = max_len > 0 ? "--maximum-length ${max_len}" : ''
    def quality = params.traditional_trim_quality ?: 0
    def adapter_args = ''

    if (params.fastp_adapter_sequence) {
        adapter_args = "--adapter '${params.fastp_adapter_sequence}'"
    } else if (adapter_fasta.name != 'NO_ADAPTER_FILE') {
        adapter_args = "--adapter file:${adapter_fasta}"
    } else if ((params.adapter_detection_method ?: 'auto') != 'auto') {
        throw new IllegalArgumentException(
            "CUTADAPT_RPF requires --fastp_adapter_sequence or --fastp_adapter_fasta"
        )
    } else {
        throw new IllegalArgumentException(
            "CUTADAPT_RPF cannot auto-detect adapters; select fastp or Trim Galore for adapter_detection_method=auto"
        )
    }

    """
    set -euo pipefail

    cutadapt \\
        --cores ${task.cpus} \\
        --minimum-length ${min_len} \\
        ${max_len_arg} \\
        --quality-cutoff ${quality} \\
        ${adapter_args} \\
        ${args} \\
        --output ${prefix}_trimmed.fastq.gz \\
        ${reads} > ${prefix}_cutadapt.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_trimmed.fastq.gz
    touch ${prefix}_cutadapt.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: 5.2
    END_VERSIONS
    """
}
