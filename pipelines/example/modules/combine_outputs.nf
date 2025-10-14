process COMBINE_OUTPUTS {
    label 'process_low'

    container "https://depot.galaxyproject.org/singularity/ubuntu:20.04"

    tag "${meta.id}"

    errorStrategy 'ignore'

    publishDir "${params.outdir}/combined",
        mode: 'copy'

    input:
    tuple val(meta), path(input_files)

    output:
    tuple val(meta), path("${meta.id}_combined.txt"),   emit: combined
    path "versions.yml",                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat ${input_files} > ${meta.id}_combined.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(cat --version | head -n1 | sed 's/cat (GNU coreutils) //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_combined.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: 8.32
    END_VERSIONS
    """
}
