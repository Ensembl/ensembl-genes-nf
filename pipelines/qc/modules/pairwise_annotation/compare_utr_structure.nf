process COMPARE_UTR_STRUCTURE {
    tag { meta.pair_id }
    label 'process_low'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'
    publishDir { "${params.outdir}/qc/pairwise_annotation/${meta.pair_id}/structural" }, mode: 'copy', pattern: 'utr_structure_comparison.tsv'

    input:
    tuple val(meta), path(source_a_models), path(source_b_models)

    output:
    tuple val(meta), path('utr_structure_comparison.tsv'), emit: results
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    pairwise_qc.py utr --pair-id '${meta.pair_id}' --source-a '${source_a_models}' --source-b '${source_b_models}' --output utr_structure_comparison.tsv
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    printf 'pair_id\\tlocus_id\\tsource_a_transcript_id\\tsource_b_transcript_id\\tsource_a_utr_available\\tsource_b_utr_available\\texact_match\\tclassification\\n' > utr_structure_comparison.tsv
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.0.0
    END_VERSIONS
    """
}
