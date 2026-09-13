process COUNT_LINES {
    label 'process_low'

    container "docker://ubuntu@sha256:c664f8f86ed5a386b0a340d981b8f81714e21a8b9c73f658c4bea56aa179d54a"

    tag "${meta.id}"

    errorStrategy 'ignore'

    publishDir "${params.outdir}/line_counts",
        mode: 'copy',
        pattern: "*_line_count.txt"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.id}_line_count.txt"), emit: counts
    path "versions.yml",                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    wc -l ${input_file} > ${meta.id}_line_count.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wc: \$(wc --version | head -n1 | sed 's/wc (GNU coreutils) //g')
    END_VERSIONS
    """

    stub:
    """
    echo "4 ${input_file}" > ${meta.id}_line_count.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wc: 8.32
    END_VERSIONS
    """
}
