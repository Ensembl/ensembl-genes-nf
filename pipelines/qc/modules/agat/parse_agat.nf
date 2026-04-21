nextflow.enable.dsl=2

process AGAT_PARSE {

    tag { meta.id }
    label 'process_light'
    publishDir "${params.outdir}/qc/agat", mode: 'copy', overwrite: true,
        pattern: "*_genebuild.csv"


    input:
        // from AGAT_RUN_STATS
        tuple val(meta), path(stats_txt)
        // path to the ensembl-genes repo containing parse_agat.py
        val ensembl_genes_repo
        // optional explicit parser path
        val agat_parser

    output:
        tuple val(meta), path("${stats_txt.simpleName}_genebuild.csv"), emit: genebuild_csv

    script:
        def out_csv    = "${stats_txt.simpleName}_genebuild.csv"
        def parser     = agat_parser ?: "${ensembl_genes_repo}/src/python/ensembl/genes/annotation-qc/parsers/parse_agat.py"
        """
        python ${parser} \\
            --input_txt ${stats_txt} \\
            --output ${out_csv}
        """
}
