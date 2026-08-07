process AGAT_PARSE {
    tag { meta.id }
    label 'qc_parser'
    publishDir "${params.outdir}/qc/agat", mode: 'copy', overwrite: true,
        pattern: "*_genebuild.csv"

    input:
        tuple val(meta), path(stats_txt)

    output:
        tuple val(meta), path("*_agat_stats_genebuild.csv"), emit: genebuild_csv
        path 'versions.yml', emit: versions

    script:
        def stem = meta.sample ?: meta.id ?: stats_txt.simpleName.replaceFirst(/_agat_stats$/, '')
        def out_csv = "${stem}_agat_stats_genebuild.csv"
        def args = task.ext.args ?: ''
        """
        annotation-qc parse-agat \
            --input_txt ${stats_txt} \
            --output ${out_csv} \
            ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            qc_parser: annotation-qc
        END_VERSIONS
        """

    stub:
        def stem = meta.sample ?: meta.id ?: stats_txt.simpleName.replaceFirst(/_agat_stats$/, '')
        """
        touch ${stem}_agat_stats_genebuild.csv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            qc_parser: annotation-qc
        END_VERSIONS
        """
}
