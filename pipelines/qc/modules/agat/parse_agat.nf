
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
        val  feature_levels_yaml

    output:
        tuple val(meta), path("${meta.sample ?: meta.id ?: stats_txt.simpleName.replaceFirst(/_agat_stats$/, '')}_agat_stats_genebuild.csv"), emit: genebuild_csv
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def stem = meta.sample ?: meta.id ?: stats_txt.simpleName.replaceFirst(/_agat_stats$/, '')
        def out_csv = "${stem}_agat_stats_genebuild.csv"
        def parser = agat_parser ?: "${ensembl_genes_repo}/src/python/ensembl/genes/annotation_qc/parsers/parse_agat.py"
        def args = task.ext.args ?: ''
        """
        python ${parser} \\
            --input_txt ${stats_txt} \\
            --output ${out_csv} \\
            ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | awk '{print \$2}')
        END_VERSIONS
        """

    stub:
        def stem = meta.sample ?: meta.id ?: stats_txt.simpleName.replaceFirst(/_agat_stats$/, '')
        """
        touch ${stem}_agat_stats_genebuild.csv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | awk '{print \$2}')
        END_VERSIONS
        """
}
