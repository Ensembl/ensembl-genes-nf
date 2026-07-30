
process PROCESS_INTERPRO {

    tag { meta.id }
    label 'process_light'
    publishDir "${params.outdir}/qc/interpro", mode: 'copy', overwrite: true,


    input:
        // from RUN_INTERPRO
        tuple val(meta), path(interpro_tsv)
        // path to the ensembl-genes repo containing parse_agat.py
        val ensembl_genes_repo

    output:
        tuple val(meta), path("${meta}_interpro_stats.tsv"), emit: interpro_stats_tsv

    script:
        def out_csv    = "${meta}_interpro_stats.tsv"
        def parser     = interpro_parser ?: "${ensembl_genes_repo}/src/python/ensembl/genes/annotation-qc/parsers/interpro.py"
        """
        python ${parser} \\
            --input_txt ${interpro_tsv} \\
            --output ${out_csv}
        """
}
