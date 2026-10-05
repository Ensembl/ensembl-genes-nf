process VALIDATE_BACKEND_SHARD_PLAN {
    tag "${backend}:${meta.id}:shard-status"
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), val(backend), path(plan), val(skip_keys)
    val success_keys

    output:
    tuple val(meta), val(backend), path("${backend}_${meta.id}_backend_shard_status.tsv"), emit: status
    path 'versions.yml', emit: versions

    script:
    """
    validate_backend_shard_plan.py ${plan} '${backend}' '${meta.id}' '${success_keys}' \\
        ${backend}_${meta.id}_backend_shard_status.tsv
    printf '"%s":\\n    backend_shard_status: python\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    validate_backend_shard_plan.py ${plan} '${backend}' '${meta.id}' '${success_keys}' \\
        ${backend}_${meta.id}_backend_shard_status.tsv
    printf '"%s":\\n    backend_shard_status: stub\\n' '${task.process}' > versions.yml
    """
}
