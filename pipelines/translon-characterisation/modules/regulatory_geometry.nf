process REGULATORY_GEOMETRY {
    label 'process_low'
    container 'oras://community.wave.seqera.io/library/python:3.12-slim'
    tag "${meta.id}"
    input:
    tuple val(meta), path(instances)
    output:
    tuple val(meta), path('*.regulatory_geometry.axis.jsonl'), emit: axis
    path 'versions.yml', emit: versions
    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    regulatory_geometry.py --instances ${instances} --output ${prefix}.regulatory_geometry.axis.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        regulatory_geometry: 0.1.0
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.regulatory_geometry.axis.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        regulatory_geometry: 0.1.0
    END_VERSIONS
    """
}
