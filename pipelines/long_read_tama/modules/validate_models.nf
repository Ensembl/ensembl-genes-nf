process VALIDATE_LONG_READ_MODELS {
    publishDir path: { meta.id == params.cohort_id ? "${params.outdir}/per_cohort/${params.cohort_id}/qc" : "${params.outdir}/per_run/${meta.id}/qc" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${backend}:${meta.id}:models"
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), val(backend), path(model_bed)

    output:
    tuple val(meta), val(backend), path("${backend}_${meta.id}_model_validation.tsv"), emit: report
    path 'versions.yml', emit: versions

    script:
    """
    validate_tama_bed.py ${model_bed} ${backend}_${meta.id}_model_validation.tsv
    printf '"%s":\\n    validator: python\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'models\\t1\\n' > ${backend}_${meta.id}_model_validation.tsv
    printf '"%s":\\n    validator: stub\\n' '${task.process}' > versions.yml
    """
}
