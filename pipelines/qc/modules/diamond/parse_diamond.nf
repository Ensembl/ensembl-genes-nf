process DIAMOND_PARSE {
    tag { meta.id }
    label 'qc_parser'
    publishDir "${params.outdir}/qc/diamond", mode: 'copy', overwrite: true

    input:
        tuple val(meta), path(genome_fasta), path(gff3), path(diamond_hits)

    output:
        tuple val(meta), path("*_translation_validity.csv"), emit: stats
        path 'versions.yml', emit: versions

    script:
        def out_dir = "${meta.id}_translation_validity"
        def out_csv = "${meta.id}_translation_validity.csv"
        """
        annotation-qc translation-validity \
            --annotation ${gff3} \
            --genome ${genome_fasta} \
            --diamond-hits ${diamond_hits} \
            --outdir ${out_dir}

        mv ${out_dir}/translation_validity_results.csv ${out_csv}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            qc_parser: \$(python -c 'import importlib.metadata; print(importlib.metadata.version("ensembl-genes"))')
        END_VERSIONS
        """

    stub:
        """
        touch ${meta.id}_translation_validity.csv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            qc_parser: \$(python -c 'import importlib.metadata; print(importlib.metadata.version("ensembl-genes"))')
        END_VERSIONS
        """
}
