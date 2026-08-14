process LOAD_REFSEQ {
    label 'process_high_memory'

    tag "${meta.id}"
    publishDir { "${params.outdir}/refseq/${meta.id}" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta),
          path(gff3),
          path(fasta),
          path(asm_report)

    tuple val(db_host),
          val(db_port),
          val(db_user),
          val(db_password)

    output:
    tuple val(meta), path("${meta.id}.load_refseq.done"), emit: loaded
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def script = "${params.ensembl_genes_repo}/src/python/ensembl/genes/ensembl_loading/gff_cli.py"


        """
        python ${script} \
            --log-file gff-loader.log \
            create-core \
            ${gff3} \
            ${fasta} \
            ${asm_report} \
            ${args} \
            --assembly-acc ${meta.id} \
            --species-name "${meta.species}" \
            --db-host ${db_host} \
            --db-port ${db_port} \
            --db-user ${db_user} \
            --db-password '${db_password}' \
            --source refseq

        touch ${meta.id}.load_refseq.done

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | awk '{print \$2}')
            gff_cli.py: unknown
        END_VERSIONS
        """

    stub:
        """
        touch ${meta.id}.load_refseq.done

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: stub
            gff_cli.py: stub
        END_VERSIONS
        """


}
