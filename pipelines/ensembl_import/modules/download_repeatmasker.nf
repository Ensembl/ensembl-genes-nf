process DOWNLOAD_REPEATMASKER {
    label 'process_medium'

    tag "${meta.id}"
    publishDir { "${params.outdir}/refseq/${meta.id}/repeatmasker" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(refseq_loaded)

    output:
    tuple val(meta), path("repeatmasker_data/**/*_rm.out.gz"), emit: repeatmasker
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    gff-loader \\
        --log-file gff-loader.log \\
        refseq \\
        download-repeatmasker \\
        --base-dir repeatmasker_data \\
        --assembly-acc ${meta.id} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print \$2}')
        ensembl-genes: \$(python -c 'from importlib.metadata import version; print(version("ensembl-genes"))')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p repeatmasker_data/${meta.id}
    touch repeatmasker_data/${meta.id}/${meta.id}_rm.out.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
        ensembl-genes: stub
    END_VERSIONS
    """
}
