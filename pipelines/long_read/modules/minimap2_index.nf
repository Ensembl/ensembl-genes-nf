process MINIMAP2_INDEX {
    tag "${meta.id}"
    label 'process_high_memory'

    conda "bioconda::minimap2=2.28"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.28--he4a0461_3' :
        'biocontainers/minimap2:2.28--he4a0461_3' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.mmi"), emit: index
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    minimap2 \\
        -d ${prefix}.mmi \\
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.mmi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: 2.28
    END_VERSIONS
    """
}
