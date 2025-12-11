process TOOL_C {
    label 'process_low'

    container "https://depot.galaxyproject.org/singularity/ubuntu:20.04"

    tag "${meta.id}"

    errorStrategy 'ignore'

    publishDir "${params.outdir}/tool_c",
        mode: 'copy'

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.id}_C.txt"),          emit: results
    path "versions.yml",                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "C output" > ${meta.id}_C.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool_c: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_C.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool_c: 1.0.0
    END_VERSIONS
    """
}
