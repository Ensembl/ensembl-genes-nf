process SUMMARISE_MODEL_COUNTS {
    tag { meta.pair_id }
    label 'process_low'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'
    publishDir { "${params.outdir}/qc/pairwise_annotation/${meta.pair_id}/structural" }, mode: 'copy', pattern: 'model_counts.tsv'

    input:
    tuple val(meta), path(source_a_models), path(source_b_models)

    output:
    tuple val(meta), path('model_counts.tsv'), emit: counts
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    pairwise_qc.py counts --pair-id '${meta.pair_id}' --source-a '${source_a_models}' --source-b '${source_b_models}' --output model_counts.tsv
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.1.0
    END_VERSIONS
    """

    stub:
    """
    printf 'pair_id\\tscope\\tsummary_level\\tlocus_count\\ttranscript_count\\texon_count\\tdistinct_locus_count\\tstatus\\n' > model_counts.tsv
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.1.0
    END_VERSIONS
    """
}
