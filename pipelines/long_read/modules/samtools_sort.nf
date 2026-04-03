process SAMTOOLS_SORT {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::samtools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0' :
        'biocontainers/samtools:1.20--h50ea8bc_0' }"

    input:
    tuple val(meta),  path(bam)
    tuple val(meta2), path(fasta), path(fai)
    val   index_format

    output:
    tuple val(meta), path("*.bam"),            optional: true, emit: bam
    tuple val(meta), path("*.cram"),           optional: true, emit: cram
    tuple val(meta), path("*.{bai,csi,crai}"), optional: true, emit: index
    path "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix   = task.ext.prefix ?: meta.id
    def args     = task.ext.args   ?: ''
    def ref_flag = fasta           ? "--reference ${fasta}" : ''
    def idx_flag = index_format    ? "--write-index" : ''
    def out_ext  = fasta           ? (index_format == 'crai' ? 'cram' : 'bam') : 'bam'
    """
    samtools sort \\
        ${args} \\
        -@ ${task.cpus} \\
        ${ref_flag} \\
        ${idx_flag} \\
        -o ${prefix}.sorted.${out_ext} \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.sorted.bam ${prefix}.sorted.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.20
    END_VERSIONS
    """
}
