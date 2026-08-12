process PAIRWISE_COUNT_GFF_TRANSCRIPTS {
    tag { "${meta.id}_${source_label}" }
    label 'process_low'

    container "https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1"

    publishDir "${params.outdir}/qc/pairwise_annotation",
        mode: 'copy',
        pattern: "*_gene_transcript_counts.tsv"

    input:
        tuple val(meta), val(source_label), path(gff)

    output:
        tuple val(meta), val(source_label), path("${meta.id}_${source_label}_gene_transcript_counts.tsv"), emit: counts
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        set -euo pipefail

        count_gff_transcripts.py \\
            --gff ${gff} \\
            --output ${meta.id}_${source_label}_gene_transcript_counts.tsv \\
            --assembly-accession ${meta.id}

        normalise_pairwise_tsv.py --input ${meta.id}_${source_label}_gene_transcript_counts.tsv --output ${meta.id}_${source_label}_gene_transcript_counts.normalised.tsv
        mv ${meta.id}_${source_label}_gene_transcript_counts.normalised.tsv ${meta.id}_${source_label}_gene_transcript_counts.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """

    stub:
        """
        printf "assembly_accession\\tsample_name\\tgene_id\\tbiotype\\tn_transcripts\\n" > ${meta.id}_${source_label}_gene_transcript_counts.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: 3.11.0
        END_VERSIONS
        """
}
