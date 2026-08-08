process EXTRACT_RPFS {
    tag "${meta.id}"
    label 'process_high'

    conda "conda-forge::python=3.10 conda-forge::biopython"
    container "ghcr.io/jackcurragh/get-rpf:0.3.0"

    input:
    tuple val(meta), path(input_file)
    path star_index

    output:
    tuple val(meta), path("*.collapsed.fa"), emit: trimmed_collapsed
    tuple val(meta), path("*.seqspec.yaml"), emit: seqspec
    tuple val(meta), path("*.extraction_report.json"), emit: report
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Use the recommended 'extract' command (alignment-based extraction)
    # Optimized with --collapsed-only for speed and disk space
    getRPF extract \\
        ${input_file} \\
        ${prefix}_trimmed.fastq \\
        -f fastq \\
        --generate-seqspec \\
        --output-format json \\
        --star-index ${star_index} \\
        --star-threads ${task.cpus} \\
        --collapsed-only \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getRPF: \$(getRPF --version 2>&1 | head -n1 | sed 's/^.*version //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_trimmed.collapsed.fa
    touch ${prefix}_trimmed.seqspec.yaml
    touch ${prefix}_trimmed.extraction_report.json
    touch ${prefix}_trimmed.report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getRPF: 0.3.0
    END_VERSIONS
    """
}
