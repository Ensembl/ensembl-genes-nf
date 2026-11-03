process COMPARE_TRANSCRIPT_FEATURES {
    tag { meta.pair_id }
    label 'process_low'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'
    publishDir { "${params.outdir}/qc/pairwise_transcripts/${meta.pair_id}/comparison" }, mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(meta), path(source_a_models), path(source_b_models), path(evidence_a), path(evidence_b)

    output:
    tuple val(meta), path('model_presence.tsv'), path('feature_presence.tsv'), path('terminal_exon_comparison.tsv'), path('structural_summary.tsv'), emit: results

    script:
    """
    set -euo pipefail
    pairwise_qc.py compare-features --source-a '${source_a_models}' --source-b '${source_b_models}' --evidence-a '${evidence_a}' --evidence-b '${evidence_b}' --pair-id '${meta.pair_id}' --model-presence model_presence.tsv --feature-presence feature_presence.tsv --terminal-exons terminal_exon_comparison.tsv --summary structural_summary.tsv
    """

    stub:
    """
    pairwise_qc.py compare-features --source-a '${source_a_models}' --source-b '${source_b_models}' --evidence-a '${evidence_a}' --evidence-b '${evidence_b}' --pair-id '${meta.pair_id}' --model-presence model_presence.tsv --feature-presence feature_presence.tsv --terminal-exons terminal_exon_comparison.tsv --summary structural_summary.tsv
    """
}
