process PAIRWISE_GFF_FEATURE_METRICS {
    tag { meta.id }
    label 'process_low'

    container "python:3.11-slim"

    publishDir "${params.outdir}/qc/pairwise_annotation/${meta.id}/feature_metrics",
        mode: 'copy',
        pattern: "*.tsv"

    input:
        tuple val(meta), path(source_a_gff), path(source_b_gff)

    output:
        tuple val(meta), path("feature_counts.tsv"), emit: features
        tuple val(meta), path("gene_metrics.tsv"), emit: gene_metrics
        tuple val(meta), path("tx_metrics.tsv"), emit: tx_metrics
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        set -euo pipefail

        compare_gff_features.py \\
            --ensembl-gff ${source_a_gff} \\
            --cat-gff ${source_b_gff} \\
            --output-dir . \\
            --assembly-accession ${meta.id} \\
            --sample-name ${meta.id}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """

    stub:
        """
        printf "assembly_accession\\tsample_name\\tsource\\tn_genes\\tn_transcripts\\tn_exon\\tn_cds\\tn_utr5\\tn_utr3\\n" > feature_counts.tsv
        printf "assembly_accession\\tsample_name\\tsource\\tgene_id\\tbiotype\\tbiotype_group\\tn_transcripts\\tgene_span_bp\\n" > gene_metrics.tsv
        printf "assembly_accession\\tsample_name\\tsource\\tgene_id\\ttranscript_id\\tbiotype\\tbiotype_group\\tn_exons\\ttx_len_bp\\tcds_len_bp\\n" > tx_metrics.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: 3.11.0
        END_VERSIONS
        """
}
