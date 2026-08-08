process FANBACK_MSA {
    label 'process_light'
    container 'oras://community.wave.seqera.io/library/python:3.12-slim'
    tag "${meta.id}"
    input:
    tuple val(meta), path(instances)
    tuple val(msa_meta), path(msa_records)
    output:
    tuple val(meta), path('*.msa.fanback.jsonl'), emit: axis
    path 'versions.yml', emit: versions
    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    fanback_axis.py --instances ${instances} --axis-records ${msa_records} --output ${prefix}.msa.fanback.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fanback_axis: 0.1.0
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.msa.fanback.jsonl
    printf '"${task.process}":\\n  fanback_axis: 0.1.0\\n' > versions.yml
    """
}
