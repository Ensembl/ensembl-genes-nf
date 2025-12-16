process EXTRACT_RPFS {
    tag "${meta.id}"
    label 'process_medium'

    conda "conda-forge::python=3.10 conda-forge::biopython"
    container "ghcr.io/jackcurragh/get-rpf:latest"

    publishDir "${params.outdir}/getRPF/extract", mode: 'copy', pattern: "*.{seqspec.yaml,extraction_report.json,report.html}"

    input:
    tuple val(meta), path(input_file)
    path star_index

    output:
    tuple val(meta), path("*_rpfs.fastq"), emit: rpfs
    tuple val(meta), path("*.seqspec.yaml"), emit: seqspec
    tuple val(meta), path("*.extraction_report.json"), emit: report
    tuple val(meta), path("*.report.html"), emit: html_report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Ensure format is specified or inferred. Usually input_file identifies it.
    // Assuming fastq input from pipeline.
    """
    getRPF extract-rpf \\
        ${input_file} \\
        ${prefix}_rpfs.fastq \\
        -f fastq \\
        --generate-seqspec \\
        --output-format json \\
        --star-index ${star_index} \\
        --star-threads ${task.cpus} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getRPF: \$(getRPF --version 2>&1 | head -n1 | sed 's/^.*version //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_rpfs.fastq
    touch ${prefix}_rpfs.seqspec.yaml
    touch ${prefix}_rpfs.extraction_report.json
    touch ${prefix}_rpfs.report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getRPF: 0.1.0
    END_VERSIONS
    """
}
