process GET_METADATA {
    label 'process_low'

    tag "${meta.id}"

    publishDir "${params.outdir}/metadata/${meta.id}", mode: 'copy', overwrite: true

    input:
    val meta
    tuple val(db_host),
          val(db_port),
          val(db_user),
          val(db_password)

    output:
    tuple val(meta),
        path("metadata/*.sql"),
        emit: sql
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def script = "${params.ensembl_genes_repo}/src/python/ensembl/genes/metadata/core_meta_data.py"

    """
    python ${script} \
        --output_dir metadata \
        --db_name ${meta.db_name} \
        --host ${db_host} \
        --port ${db_port} \
        --team "genebuild" \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print \$2}')
        core_metadata.py: unknown
    END_VERSIONS
    """

    stub:
    """
    mkdir -p metadata
    touch metadata/${meta.db_name}.sql

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
        core_metadata.py: stub
    END_VERSIONS
    """
}
