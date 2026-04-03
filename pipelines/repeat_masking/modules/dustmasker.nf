// DUSTMASKER — NCBI low-complexity sequence masker
// Outputs a BED file of low-complexity regions.

process DUSTMASKER {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::blast=2.15"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.15.0--pl5321h6f7f691_1' :
        'biocontainers/blast:2.15.0--pl5321h6f7f691_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.dust.bed"), emit: bed
    path "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args   = task.ext.args   ?: ''
    """
    dustmasker \\
        -in    ${fasta} \\
        -outfmt interval \\
        ${args} \\
    | dust_to_bed.py --seqid-from-fasta ${fasta} > ${prefix}.dust.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dustmasker: \$(dustmasker -version 2>&1 | head -1 || echo 'unknown')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.dust.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dustmasker: 2.15.0
    END_VERSIONS
    """
}
