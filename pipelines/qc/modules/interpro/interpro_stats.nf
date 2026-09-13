process PROCESS_INTERPRO {
    tag { meta.id }
    label 'qc_parser'
    publishDir "${params.outdir}/qc/interpro", mode: 'copy', overwrite: true

    input:
        tuple val(meta), path(interpro_tsv), path(query_protein)

    output:
        tuple val(meta), path("*_hits.tsv"), emit: interpro_stats_tsv
        path 'versions.yml', emit: versions

    script:
        def stem = meta.sample ?: meta.id ?: interpro_tsv.simpleName.replaceFirst(/\.tsv$/, '')
        def out_dir = "${stem}_interpro_stats"
        def out_tsv = "${stem}_hits.tsv"
        """
        annotation-qc parse-interpro \
            --input_tsv ${interpro_tsv} \
            --output ${out_dir} \
            --query_protein ${query_protein}

        mv ${out_dir}/interpro_hit_summary.tsv ${out_tsv}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            qc_parser: \$(python -c 'import importlib.metadata; print(importlib.metadata.version("ensembl-genes"))')
        END_VERSIONS
        """

    stub:
        def stem = meta.sample ?: meta.id ?: interpro_tsv.simpleName.replaceFirst(/\.tsv$/, '')
        """
        touch ${stem}_hits.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            qc_parser: \$(python -c 'import importlib.metadata; print(importlib.metadata.version("ensembl-genes"))')
        END_VERSIONS
        """
}
