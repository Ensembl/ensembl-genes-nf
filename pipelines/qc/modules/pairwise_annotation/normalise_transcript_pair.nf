process NORMALISE_TRANSCRIPT_PAIR {
    tag { meta.pair_id }
    label 'process_low'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'
    publishDir { "${params.outdir}/qc/pairwise_transcripts/${meta.pair_id}/models" }, mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(meta), path(source_a_annotation), path(source_b_annotation)

    output:
    tuple val(meta), path('source_a.models.tsv'), path('source_b.models.tsv'), path('model_conversion_audit.tsv'), emit: models

    script:
    """
    set -euo pipefail
    pairwise_qc.py normalise-pair --source-a '${source_a_annotation}' --source-b '${source_b_annotation}' --format-a '${meta.format_a}' --format-b '${meta.format_b}' --pair-id '${meta.pair_id}' --assembly '${meta.assembly}' --output-a source_a.models.tsv --output-b source_b.models.tsv --audit model_conversion_audit.tsv
    """

    stub:
    """
    pairwise_qc.py normalise-pair --source-a '${source_a_annotation}' --source-b '${source_b_annotation}' --format-a '${meta.format_a}' --format-b '${meta.format_b}' --pair-id '${meta.pair_id}' --assembly '${meta.assembly}' --output-a source_a.models.tsv --output-b source_b.models.tsv --audit model_conversion_audit.tsv
    """
}
