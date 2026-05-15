process PAIRWISE_TRANSCRIPT_CONCORDANCE {
    tag { meta.id }
    label 'process_low'

    container "python:3.11-slim"

    publishDir "${params.outdir}/qc/pairwise_annotation/${meta.id}/transcript_concordance",
        mode: 'copy',
        pattern: "*_transcript_concordance.tsv"

    input:
        tuple val(meta), path(source_a_gff), path(source_b_gff), path(rbh_pairs)

    output:
        tuple val(meta), path("${meta.id}_transcript_concordance.tsv"), emit: metrics
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        set -euo pipefail

        calculate_transcript_concordance.py \\
            --ensembl-gff ${source_a_gff} \\
            --cat-gff ${source_b_gff} \\
            --pairs ${rbh_pairs} \\
            --output ${meta.id}_transcript_concordance.tsv \\
            --assembly-accession ${meta.id} \\
            --sample-name ${meta.id}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """

    stub:
        """
        printf "assembly_accession\\tsample_name\\tensembl_gene_id\\tcat_gene_id\\tensembl_biotype\\tn_ensembl_transcripts\\tn_cat_transcripts\\ttranscript_concordance_rate\\n" > ${meta.id}_transcript_concordance.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: 3.11.0
        END_VERSIONS
        """
}
