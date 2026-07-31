process LOAD_REFSEQ {
    tag "${meta.id}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(meta),
          path(gff3),
          path(fasta),
          path(asm_report)

    tuple val(db_host),
          val(db_port),
          val(db_user),
          val(db_password)

    script:
        def script = "${ensembl_genes_repo}/src/python/ensembl/genes/ensembl_loading/gff_cli.py"


        """
        python ${script} \
            gff-loader \
            --log-file gff-loader.log \
            create-core \
            ${gff3} \
            ${fasta} \
            ${asm_report} \
            --assembly-acc ${meta.id} \
            --species-name ${meta.species} \
            --db-host ${db_host} \
            --db-port ${db_port} \
            --db-user ${db_user} \
            --db-password ${db_password} \
            --source refseq
        """



}
