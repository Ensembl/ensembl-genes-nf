process NORMALISE_TRANSCRIPT_MODELS {
    tag { meta.pair_id }
    label 'process_low'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'
    publishDir { "${params.outdir}/qc/pairwise_annotation/${meta.pair_id}/models" }, mode: 'copy', pattern: '*.models.tsv'

    input:
    tuple val(meta), path(source_a_annotation), path(source_b_annotation)

    output:
    tuple val(meta), path('source_a.models.tsv'), path('source_b.models.tsv'), emit: models
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    pairwise_qc.py normalise --input '${source_a_annotation}' --source 'source_a' \\
      --pair-id '${meta.pair_id}' --assembly '${meta.assembly}' --format '${meta.format_a}' \\
      --output source_a.models.tsv
    pairwise_qc.py normalise --input '${source_b_annotation}' --source 'source_b' \\
      --pair-id '${meta.pair_id}' --assembly '${meta.assembly}' --format '${meta.format_b}' \\
      --output source_b.models.tsv
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    printf 'pair_id\\tsource\\toriginal_gene_id\\toriginal_transcript_id\\tcontig\\tstrand\\ttranscript_start\\ttranscript_end\\texon_count\\texon_coordinates\\tintron_chain\\tfive_prime_end\\tthree_prime_end\\tassembly\\tsource_format\\tmodel_status\\tutr_coordinates\\tcds_coordinates\\n' > source_a.models.tsv
    cp source_a.models.tsv source_b.models.tsv
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.0.0
    END_VERSIONS
    """
}
