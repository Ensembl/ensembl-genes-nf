process FASTP {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::fastp=0.23.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--h5f740d0_0' :
        'biocontainers/fastp:0.23.4--h5f740d0_0' }"

    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: '*.json'
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: '*.html'
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: '*_clipped.fastq.gz'

    input:
    tuple val(meta), path(reads), path(adapter_fasta)

    output:
    tuple val(meta), path("*_clipped_provided.fastq.gz"), emit: trimmed_fastq_provided
    tuple val(meta), path("*_clipped_final.fastq.gz"), emit: trimmed_fastq
    tuple val(meta), path("*_provided_fastp.json"), emit: json_provided
    tuple val(meta), path("*_auto_fastp.json"), emit: json_auto
    tuple val(meta), path("*_provided_fastp.html"), emit: html_provided
    tuple val(meta), path("*_auto_fastp.html"), emit: html_auto
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Run with provided adapters
    fastp \\
        -i $reads \\
        -o ${prefix}_clipped_provided.fastq.gz \\
        $args \\
        --length_required 20 \\
        --adapter_fasta $adapter_fasta \\
        --json ${prefix}_provided_fastp.json \\
        --html ${prefix}_provided_fastp.html \\
        --thread $task.cpus

    # Run with automatic adapter detection on the output of the first run
    fastp \\
        -i ${prefix}_clipped_provided.fastq.gz \\
        -o ${prefix}_clipped_final.fastq.gz \\
        $args \\
        --length_required 20 \\
        --json ${prefix}_auto_fastp.json \\
        --html ${prefix}_auto_fastp.html \\
        --thread $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_clipped_provided.fastq.gz
    touch ${prefix}_clipped_final.fastq.gz
    touch ${prefix}_provided_fastp.json
    touch ${prefix}_auto_fastp.json
    touch ${prefix}_provided_fastp.html
    touch ${prefix}_auto_fastp.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: 0.23.4
    END_VERSIONS
    """
}
