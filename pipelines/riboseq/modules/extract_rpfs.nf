process EXTRACT_RPFS {
    tag "${meta.id}"
    label 'process_medium'

    conda "conda-forge::python=3.10 conda-forge::biopython"
    container "ghcr.io/jackcurragh/get-rpf:main"

    publishDir "${params.outdir}/getRPF/extract", mode: 'copy', pattern: "*.{extraction_report.json}"

    input:
    tuple val(meta), path(input_file)
    path star_index

    output:
    tuple val(meta), path("*_rpfs.fastq"), emit: rpfs
    tuple val(meta), path("*.extraction_report.json"), emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_size = params.getrpf_max_reads ?: 10000
    def preserve_umi = params.getrpf_preserve_umi ? '--preserve-umi' : ''
    """
    gzip -d -c ${input_file} > temp_input.fastq
    # Use the recommended 'extract' command (alignment-based extraction)
    getRPF extract \\
        temp_input.fastq \\
        ${prefix}_rpfs.fastq \\
        --star-index ${star_index} \\
        --sample-size ${sample_size} \\
        --threads ${task.cpus} \\
        ${preserve_umi} \\
        --output-report ${prefix}.extraction_report.json \\
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
    touch ${prefix}.extraction_report.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getRPF: 0.1.0
    END_VERSIONS
    """
}
