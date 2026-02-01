/*
 * MAKE_TRANSCRIPTOME
 * Extract transcriptome sequences from genome using GTF annotation
 */

process MAKE_TRANSCRIPTOME {
    tag "${organism}_${version}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hd03093a_1' :
        'quay.io/biocontainers/gffread:0.12.7--hd03093a_1' }"

    publishDir "${params.outdir}/organism_setup/${organism}/${version}", mode: 'copy'

    input:
    path(gtf)
    path(fasta)
    val(organism)
    val(version)

    output:
    path "*.transcripts.fa",   emit: transcripts
    path "versions.yml",       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${organism}_${version}"
    """
    gffread \\
        -w ${prefix}.transcripts.fa \\
        -g ${fasta} \\
        ${gtf} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${organism}_${version}"
    """
    touch ${prefix}.transcripts.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: 0.12.7
    END_VERSIONS
    """
}
