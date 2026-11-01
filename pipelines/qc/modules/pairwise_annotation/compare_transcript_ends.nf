process COMPARE_TRANSCRIPT_ENDS {
    tag { meta.pair_id }
    label 'process_low'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'
    publishDir { "${params.outdir}/qc/pairwise_annotation/${meta.pair_id}/structural" }, mode: 'copy', pattern: 'transcript_end_comparison.tsv'

    input:
    tuple val(meta), path(source_a_models), path(source_b_models)

    output:
    tuple val(meta), path('transcript_end_comparison.tsv'), emit: results
    path 'versions.yml', emit: versions

    script:
    def threshold = params.pairwise_end_threshold ?: 50
    """
    set -euo pipefail
    pairwise_qc.py ends --pair-id '${meta.pair_id}' --source-a '${source_a_models}' --source-b '${source_b_models}' --threshold '${threshold}' --output transcript_end_comparison.tsv
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    printf 'pair_id\\tlocus_id\\tend\\tsource_a_transcript_id\\tsource_b_locus_id\\tsource_b_transcript_id\\tdistance_bp\\tclassification\\n' > transcript_end_comparison.tsv
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.0.0
    END_VERSIONS
    """
}
