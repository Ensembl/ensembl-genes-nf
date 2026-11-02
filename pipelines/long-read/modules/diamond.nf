process PREPARE_DIAMOND_QUERY {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${meta.id}:${backend}:${scope}"
    input:
    tuple val(meta), val(backend), val(scope), path(orfs)
    output:
    tuple val(meta), val(backend), val(scope), path('diamond_query.fa'), path('diamond_query.status.tsv'), emit: query
    script:
    """
    python3 ${projectDir}/bin/prepare_diamond_query.py --input ${orfs} --output diamond_query.fa --status diamond_query.status.tsv
    """
    stub:
    """
    cp ${orfs} diamond_query.fa
    printf 'status\\tstate\\nquery\\tVALID\\n' > diamond_query.status.tsv
    """
}

process EMIT_EMPTY_DIAMOND_HITS {
    label 'process_low'
    container params.use_containers ? params.tool_containers.python_helper : null
    tag "${meta.id}:${backend}:${scope}"
    input:
    tuple val(meta), val(backend), val(scope), path(query), path(status)
    output:
    tuple val(meta), val(backend), val(scope), path('diamond.tsv'), emit: hits
    script:
    """
    python3 ${projectDir}/bin/write_text.py --base64 '' --output diamond.tsv
    """
    stub:
    """
    touch diamond.tsv
    """
}
