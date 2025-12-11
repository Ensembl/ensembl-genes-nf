process TOOL_A {
    label 'process_low'

    container "https://depot.galaxyproject.org/singularity/ubuntu:20.04"

    tag "${meta.id}"

    errorStrategy 'ignore'

    publishDir "${params.outdir}/tool_a",
        mode: 'copy'

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.id}_A.txt"),          emit: results
    path "versions.yml",                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "A output" > ${meta.id}_A.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool_a: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_A.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool_a: 1.0.0
    END_VERSIONS
    """
}
