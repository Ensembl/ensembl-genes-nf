process EXTRACT_RPFS_OPTIMIZED {
    tag "${meta.id}"
    label 'process_high'

    conda "conda-forge::python=3.10 conda-forge::biopython"
    container "ghcr.io/jackcurragh/get-rpf:main"

    publishDir "${params.outdir}/getRPF/extract", mode: 'copy', pattern: "*.extraction_report.json"
    publishDir "${params.outdir}/collapsed_fa", mode: 'copy', pattern: "*.collapsed.fa"

    input:
    tuple val(meta), path(input_file)
    path star_index

    output:
    tuple val(meta), path("*.collapsed.fa"), emit: collapsed_fasta
    tuple val(meta), path("*_rpfs.fastq"), emit: rpfs, optional: true
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
    # Use the recommended 'extract' command (alignment-based extraction)
    # Optimized with --collapsed-only for speed and disk space
    getRPF extract \\
        ${input_file} \\
        ${prefix}_rpfs.fastq \\
        --star-index ${star_index} \\
        --sample-size ${sample_size} \\
        --threads ${task.cpus} \\
        --collapsed-only \\
        ${preserve_umi} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getRPF: \$(getRPF --version 2>&1 | head -n1 | sed 's/^.*version //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_rpfs.collapsed.fa
    touch ${prefix}_rpfs.extraction_report.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getRPF: 0.1.0
    END_VERSIONS
    """
}
