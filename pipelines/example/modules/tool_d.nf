process TOOL_D {
    label 'process_low'

    container "docker://ubuntu@sha256:c664f8f86ed5a386b0a340d981b8f81714e21a8b9c73f658c4bea56aa179d54a"

    tag "${meta.id}"

    errorStrategy 'ignore'

    publishDir "${params.outdir}/tool_d",
        mode: 'copy',
        pattern: "*_D.txt"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.id}_D.txt"),          emit: results
    path "versions.yml",                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "D output" > ${meta.id}_D.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool_d: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_D.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool_d: 1.0.0
    END_VERSIONS
    """
}
