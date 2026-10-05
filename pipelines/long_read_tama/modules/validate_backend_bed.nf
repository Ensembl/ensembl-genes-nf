process VALIDATE_BACKEND_BED {
    tag "${backend}:${meta.id}:${shard}:bed12"
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), val(backend), val(shard), path(bed)

    output:
    tuple val(meta), val(backend), val(shard), path('validated_models.bed'), emit: bed
    tuple val(meta), val(backend), val(shard), path('backend_status.tsv'), emit: status
    path 'versions.yml', emit: versions

    script:
    """
    validate_tama_bed.py ${bed} backend_validation.tsv \\
        --validated-output validated_models.bed --run '${meta.id}' --shard '${shard}' \\
        --backend '${backend}' --status-output backend_status.tsv
    printf '"%s":\\n    backend_bed_validator: python\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'chrStub\\t0\\t4\\t${backend}.${meta.id}.${shard}.1\\t0\\t+\\t0\\t4\\t0\\t1\\t4,\\t0,\\n' > validated_models.bed
    printf 'backend\\trun_accession\\tshard\\tstatus\\tmodels\\n${backend}\\t${meta.id}\\t${shard}\\tSUCCESS\\t1\\n' > backend_status.tsv
    printf '"%s":\\n    backend_bed_validator: stub\\n' '${task.process}' > versions.yml
    """
}
