process SAMTOOLS_INDEX {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::samtools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0' :
        'biocontainers/samtools:1.20--h50ea8bc_0' }"

    input:
    tuple val(meta), path(sorted_bam)

    output:
    tuple val(meta), path("${sorted_bam}", arity: '1'), path("*.bai", arity: '1'), emit: bam_and_bai
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools index \\
        -@ ${task.cpus} \\
        ${args} \\
        ${sorted_bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${sorted_bam}.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.20
    END_VERSIONS
    """
}
