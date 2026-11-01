process ANNOTATE_MODEL_EVIDENCE {
    tag { "${meta.pair_id}_${source_label}" }
    label 'process_low'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'
    publishDir { "${params.outdir}/qc/pairwise_annotation/${meta.pair_id}/evidence" }, mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(meta), val(source_label), path(models), path(evidence)

    output:
    tuple val(meta), val(source_label), path("${source_label}.evidence.tsv"), emit: annotated
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    pairwise_qc.py evidence --models '${models}' --evidence '${evidence}' --output '${source_label}.evidence.tsv'
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    printf 'pair_id\\tsource\\ttranscript_id\\tread_count\\tjunction_support\\tfull_length_support\\tsource_tool\\n' > '${source_label}.evidence.tsv'
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.0.0
    END_VERSIONS
    """
}
