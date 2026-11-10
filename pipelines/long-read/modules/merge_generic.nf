process MERGE_BACKEND_MODELS {
    label 'process_medium'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${meta.id}:${backend}:${scope}"
    input:
    tuple val(meta), val(backend), path(gtfs), val(scope)
    output:
    tuple val(meta), val(backend), val(scope), path('merged.gtf'), emit: merged
    script:
    """
    python3 ${projectDir}/bin/merge_gtfs.py --inputs ${gtfs} --output merged.gtf --backend '${backend}' --status merge_gtf.status.tsv --scope '${scope}'
    """
    stub:
    """
    cat ${gtfs} > merged.gtf
    """
}
