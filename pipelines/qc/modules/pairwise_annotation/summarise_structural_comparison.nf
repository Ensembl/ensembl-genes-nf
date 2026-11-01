process SUMMARISE_STRUCTURAL_COMPARISON {
    tag { meta.pair_id }
    label 'process_low'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'
    publishDir { "${params.outdir}/qc/pairwise_annotation/${meta.pair_id}/structural" }, mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(meta), path(intron), path(ends), path(utr)

    output:
    tuple val(meta), path('structural_summary.tsv'), emit: summary
    tuple val(meta), path('disagreement_loci.tsv'), emit: disagreements
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    pairwise_qc.py summary --pair-id '${meta.pair_id}' --intron '${intron}' --ends '${ends}' --utr '${utr}' \\
      --output structural_summary.tsv --disagreements disagreement_loci.tsv
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    printf 'pair_id\\tn_intron_comparisons\\tn_exact_intron_chains\\tn_partial_intron_chains\\tn_end_comparisons\\tn_ends_within_threshold\\tn_utr_comparisons\\tn_exact_utr_matches\\n' > structural_summary.tsv
    printf 'pair_id\\tlocus_id\\tsource_a_transcript_id\\tsource_b_transcript_id\\tfeature\\tclassification\\n' > disagreement_loci.tsv
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.0.0
    END_VERSIONS
    """
}
