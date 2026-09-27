process TYPED_AXIS {
    label 'process_low'
    container 'oras://community.wave.seqera.io/library/python:3.12-slim'
    tag "${meta.id}:${axis_name}"

    input:
    tuple val(meta), val(axis_name), path(instances)

    output:
    tuple val(meta), val(axis_name), path('*.axis.jsonl'), emit: axis
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${axis_name}"
    """
    translon_characterise.py axis --axis ${axis_name} --input ${instances} --output ${prefix}.axis.jsonl ${args}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${axis_name}"
    """
    touch ${prefix}.axis.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
