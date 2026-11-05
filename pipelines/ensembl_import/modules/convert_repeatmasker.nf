process CONVERT_REPEATMASKER {
    label 'process_light'

    tag "${meta.id}"
    publishDir { "${params.outdir}/refseq/${meta.id}/repeatmasker" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(repeatmasker), path(assembly_report)

    output:
    tuple val(meta), path("${meta.id}.repeatmasker.gtf"), emit: gtf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    gff-loader \\
        --log-file gff-loader.log \\
        refseq \\
        convert-repeatmasker \\
        ${repeatmasker} \\
        ${assembly_report} \\
        --output ${meta.id}.repeatmasker.gtf \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print \$2}')
        ensembl-genes: \$(python -c 'from importlib.metadata import version; print(version("ensembl-genes"))')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.repeatmasker.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
        ensembl-genes: stub
    END_VERSIONS
    """
}
