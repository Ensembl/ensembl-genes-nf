process PREPARE_TRANSCRIPT_MODELS {
    tag "${meta.id}"
    label 'process_low'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'python:3.12-bookworm'
    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    output:
    tuple val(meta), path(bam), path(bai), path('models.genePred'), path('models.bed12'), emit: models
    script:
    """
    python3 ${params.translon_analysis_bin}/gtf_to_transcript_models.py --gtf ${gtf} --genepred models.genePred --bed12 models.bed12
    printf '"%s":\n    transcript_models: 1\n' '${task.process}' > versions.yml
    """
    stub:
    """
    touch models.genePred models.bed12
    printf '"stub":\n    transcript_models: stub\n' > versions.yml
    """
}
