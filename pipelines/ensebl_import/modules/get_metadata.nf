process GET_METADATA {

    tag "${meta.id}"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(meta)
    tuple val(db_host),
          val(db_port),
          val(db_user),
          val(db_password)

    output:
    tuple val(meta),
        path("metadata/*.sql"),
        emit: sql

    script:
    def script = "${ensembl_genes_repo}/src/python/ensembl/genes/metadata/core_metadata.py"

    """
    python ${script} \
        --output_dir metadata \
        --db_name ${meta.db_name} \
        --host ${db_host} \
        --port ${db_port} \
        --team "genebuild"
    """
}
