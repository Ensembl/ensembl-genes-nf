process PAIRWISE_TRANSCRIPT_CONCORDANCE {
    tag { meta.id }
    label 'process_low'

    container "https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1"

    publishDir "${params.outdir}/qc/pairwise_annotation",
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

        normalise_pairwise_tsv.py --input ${meta.id}_transcript_concordance.tsv --output ${meta.id}_transcript_concordance.normalised.tsv
        mv ${meta.id}_transcript_concordance.normalised.tsv ${meta.id}_transcript_concordance.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """

    stub:
        """
        printf "assembly_accession\\tsample_name\\tsource_a_gene_id\\tsource_b_gene_id\\tsource_a_biotype\\tn_source_a_transcripts\\tn_source_b_transcripts\\ttranscript_concordance_rate\\n" > ${meta.id}_transcript_concordance.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: 3.11.0
        END_VERSIONS
        """
}
