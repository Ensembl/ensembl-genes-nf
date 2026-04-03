process BEDTOOLS_MASKFASTA {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::bedtools=2.31"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_1' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_1' }"

    input:
    tuple val(meta),  path(bed)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.softmasked.fa"), emit: fasta
    path "versions.yml",                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args   = task.ext.args   ?: '-soft'
    """
    bedtools maskfasta \\
        -fi  ${fasta} \\
        -bed ${bed} \\
        -fo  ${prefix}.softmasked.fa \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.softmasked.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: 2.31.1
    END_VERSIONS
    """
}
