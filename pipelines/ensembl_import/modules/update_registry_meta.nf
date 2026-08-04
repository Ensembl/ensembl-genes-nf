process UPDATE_REGISTRY {

    tag "${meta.id}"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(meta), val(db_name)
    tupleval(db_host),
          val(db_port),
          val(db_user),
          val(db_password)

    output:
    tuple val(meta)
    
    emit: refseq

    script:
    def script = "${ensembl_genes_repo}/src/python/ensembl/genes/metadata/core_metadata.py"

    """
    python ${script} \
        --output_dir ${params.outdir} \
        --db_name ${db_name} \
        --host ${db_host} \
        --port ${db_port} \
        --team "genebuild" 
    """
}
