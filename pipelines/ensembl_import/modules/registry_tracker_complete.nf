process REGISTRY_TRACKER_COMPLETE {
    label 'python'

    tag "${meta.id}"

    input:
    val meta

    output:
    tuple val(meta), path("${meta.id}.registry_complete.done"), emit: completed
    path 'versions.yml', emit: versions

    script:
    """
    registry_tracker.py complete --import-status-id ${meta.registry_import_status_id}
    touch ${meta.id}.registry_complete.done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print \$2}')
        registry_tracker.py: unknown
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.registry_complete.done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
        registry_tracker.py: stub
    END_VERSIONS
    """
}
