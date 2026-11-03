process REGISTRY_TRACKER_START {
    label 'python'

    tag "${meta.id}"

    input:
    val meta

    output:
    tuple val(meta), path('import_status_id.txt'), emit: started
    path 'versions.yml', emit: versions

    script:
    def skipStartCheck = params.skip_start_check == 'true' ? '--skip-start-check' : ''

    """
    registry_tracker.py start \
        --accession ${meta.id} \
        --genebuilder '${params.genebuilder}' \
        ${skipStartCheck} \
        > import_status_id.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print \$2}')
        registry_tracker.py: unknown
    END_VERSIONS
    """

    stub:
    """
    echo stub-${meta.id} > import_status_id.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
        registry_tracker.py: stub
    END_VERSIONS
    """
}
