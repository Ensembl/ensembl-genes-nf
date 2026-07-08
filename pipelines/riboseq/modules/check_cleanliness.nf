process CHECK_CLEANLINESS {
    tag "${meta.id}"
    label 'process_medium'

    conda "conda-forge::python=3.10 conda-forge::biopython"
    container "ghcr.io/jackcurragh/get-rpf:0.2.4"

    publishDir "${params.outdir}/getRPF/check", mode: 'copy', pattern: "*_{report,rpf_checks}.txt"

    input:
    tuple val(meta), path(input_file)
    val count_pattern

    output:
    tuple val(meta), path("*_report.txt"), emit: report
    tuple val(meta), path("*rpf_checks.txt"), emit: rpf_checks
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    getRPF check-cleanliness ${input_file} \\
        --format collapsed \\
        --output . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getRPF: \$(getRPF --version 2>&1 | sed 's/getRPF version //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_report.txt
    touch ${prefix}_rpf_checks.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getRPF: 0.2.4
    END_VERSIONS
    """
}
