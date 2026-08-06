
process PROCESS_INTERPRO {

    tag { meta.id }
    label 'process_light'
    publishDir "${params.outdir}/qc/interpro", mode: 'copy', overwrite: true

    input:
        tuple val(meta), path(interpro_tsv), path(query_protein)
        val ensembl_genes_repo
        val interpro_parser

    output:
        tuple val(meta), path("${meta.sample ?: meta.id ?: interpro_tsv.simpleName.replaceFirst(/\\.tsv$/, '')}_hits.tsv"), emit: interpro_stats_tsv
        path "versions.yml", emit: versions

    script:
        def stem   = meta.sample ?: meta.id ?: interpro_tsv.simpleName.replaceFirst(/\\.tsv$/, '')
        def out_dir = "${stem}_interpro_stats"
        def out_tsv = "${stem}_hits.tsv"
        def parser = interpro_parser ?: "${ensembl_genes_repo}/src/python/ensembl/genes/annotation_qc/parsers/interpro.py"
        """
        python ${parser} \\
            --input_tsv ${interpro_tsv} \\
            --output ${out_dir} \\
            --query_protein ${query_protein}

        mv ${out_dir}/interpro_hit_summary.tsv ${out_tsv}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | awk '{print \$2}')
        END_VERSIONS
        """

    stub:
        def stem = meta.sample ?: meta.id ?: interpro_tsv.simpleName.replaceFirst(/\\.tsv$/, '')
        """
        touch ${stem}_hits.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: 3.11.0
        END_VERSIONS
        """
}
