process LOAD_REPEATMASKER {
    label 'process_high_memory'

    tag "${meta.id}"
    publishDir { "${params.outdir}/refseq/${meta.id}/repeatmasker" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(gtf)

    tuple val(db_host),
          val(db_port),
          val(db_user),
          val(db_password)

    output:
    tuple val(meta), path("${meta.id}.load_repeatmasker.done"), emit: loaded
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    gff-loader \\
        --log-file gff-loader.log \\
        load-single-line-features \\
        ${gtf} \\
        --analysis-name repeatmasker \\
        --db-name ${meta.db_name} \\
        --db-host ${db_host} \\
        --db-port ${db_port} \\
        --db-user ${db_user} \\
        --db-password '${db_password}' \\
        ${args}

    touch ${meta.id}.load_repeatmasker.done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print \$2}')
        ensembl-genes: \$(python -c 'from importlib.metadata import version; print(version("ensembl-genes"))')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.load_repeatmasker.done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
        ensembl-genes: stub
    END_VERSIONS
    """
}
