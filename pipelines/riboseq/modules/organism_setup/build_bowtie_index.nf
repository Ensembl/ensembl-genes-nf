/*
 * BUILD_BOWTIE_INDEX
 * Build Bowtie1 index optimized for short Ribo-seq reads (20-35nt)
 */

process BUILD_BOWTIE_INDEX {
    tag "${index_name}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie:1.3.1--py39hd16f23e_2' :
        'quay.io/biocontainers/bowtie:1.3.1--py39hd16f23e_2' }"

    publishDir "${params.outdir}/organism_setup/${organism}/${version}", mode: 'copy'

    input:
    path(fasta)
    val(index_name)
    val(organism)
    val(version)

    output:
    path "${index_name}_index",   emit: index
    path "versions.yml",          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Optimized for short Ribo-seq reads (20-35nt)
    // --offrate 0: maximum sensitivity for short reads (larger index)
    // --ftabchars 5: smaller initial lookup table for short seeds
    """
    mkdir -p ${index_name}_index

    bowtie-build \\
        ${args} \\
        --threads ${task.cpus} \\
        --offrate 0 \\
        --ftabchars 5 \\
        ${fasta} \\
        ${index_name}_index/${fasta.baseName}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: \$(echo \$(bowtie --version 2>&1) | sed 's/^.*bowtie-align-s version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${index_name}_index
    touch ${index_name}_index/${fasta.baseName}.1.ebwt
    touch ${index_name}_index/${fasta.baseName}.2.ebwt
    touch ${index_name}_index/${fasta.baseName}.3.ebwt
    touch ${index_name}_index/${fasta.baseName}.4.ebwt
    touch ${index_name}_index/${fasta.baseName}.rev.1.ebwt
    touch ${index_name}_index/${fasta.baseName}.rev.2.ebwt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: 1.3.1
    END_VERSIONS
    """
}
