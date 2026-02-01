/*
 * GGET_DOWNLOAD
 * Download genome and annotation files from Ensembl using gget
 */

process GGET_DOWNLOAD {
    tag "${organism}_${ensembl_version}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gget:0.29.0--pyhdfd78af_0' :
        'quay.io/biocontainers/gget:0.29.0--pyhdfd78af_0' }"

    publishDir "${params.outdir}/organism_setup/${organism}/${ensembl_version}", mode: 'copy'

    input:
    val(organism)
    val(ensembl_version)
    val(which)

    output:
    path "*.fa",           emit: genome_fasta, optional: true
    path "*.gtf",          emit: genome_gtf, optional: true
    path "versions.yml",   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def which_arg = which ? "--which ${which.join(',')}" : '--which dna,gtf'
    def release_arg = ensembl_version ? "--release ${ensembl_version}" : ''
    """
    gget ref ${organism} \\
        ${which_arg} \\
        ${release_arg} \\
        ${args} \\
        --download

    # Decompress downloaded files
    gzip -d *.gz || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gget: \$(gget --version | sed 's/gget //')
    END_VERSIONS
    """

    stub:
    """
    touch genome.fa
    touch annotation.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gget: 0.29.0
    END_VERSIONS
    """
}
