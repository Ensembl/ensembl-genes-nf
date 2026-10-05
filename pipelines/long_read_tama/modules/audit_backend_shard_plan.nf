process AUDIT_BACKEND_SHARD_PLAN {
    tag "${backend}:${meta.id}:shard-plan"
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), val(backend), path(plan), val(skip_keys)

    output:
    tuple val(meta), val(backend), path("${backend}_${meta.id}_backend_shard_plan.tsv"), emit: plan
    path 'versions.yml', emit: versions

    script:
    """
    audit_backend_shard_plan.py ${plan} '${backend}' '${meta.id}' '${skip_keys}' ${backend}_${meta.id}_backend_shard_plan.tsv
    printf '"%s":\\n    backend_shard_plan_audit: python\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    audit_backend_shard_plan.py ${plan} '${backend}' '${meta.id}' '${skip_keys}' ${backend}_${meta.id}_backend_shard_plan.tsv
    printf '"%s":\\n    backend_shard_plan_audit: stub\\n' '${task.process}' > versions.yml
    """
}
