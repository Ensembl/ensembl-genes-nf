process EVOLUTIONARY_HISTORY {
    label 'process_light'
    container 'oras://community.wave.seqera.io/library/python:3.12-slim'
    tag "${meta.id}"

    input:
    tuple val(meta), path(instances)
    path evidence

    output:
    tuple val(meta), path('*.evolutionary_history.axis.jsonl'), emit: axis
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    evolutionary_history.py --instances ${instances} --evidence ${evidence} \\
        --output ${prefix}.evolutionary_history.axis.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        evolutionary_history: read-only-candidate-evidence-v1
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    printf '{"instance_id":"stub","axis":{"evolutionary_history":{"state":"uninformative","reason":"stub"}}}\n' > ${prefix}.evolutionary_history.axis.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
        evolutionary_history: read-only-candidate-evidence-v1
    END_VERSIONS
    """
}
