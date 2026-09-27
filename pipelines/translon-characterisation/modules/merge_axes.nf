process MERGE_AXES {
    label 'process_low'
    container 'oras://community.wave.seqera.io/library/python:3.12-slim'
    tag "${meta.id}"

    input:
    tuple val(meta), path(axis_files)

    output:
    tuple val(meta), path('*.characterised.jsonl'), emit: instances
    path 'versions.yml', emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    translon_characterise.py merge-axes --inputs ${axis_files} --output ${prefix}.characterised.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.characterised.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
