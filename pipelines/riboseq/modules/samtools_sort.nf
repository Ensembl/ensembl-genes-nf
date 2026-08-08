process SAMTOOLS_SORT {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::samtools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0' :
        'biocontainers/samtools:1.20--h50ea8bc_0' }"

    // The sorted BAM is immediately indexed and consumed downstream. Keep it
    // in the work graph rather than publishing a second intermediate copy.

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.sorted.bam"), emit: bam
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools sort \\
        -@ ${task.cpus} \\
        -o ${prefix}.sorted.bam \\
        ${args} \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.20
    END_VERSIONS
    """
}
