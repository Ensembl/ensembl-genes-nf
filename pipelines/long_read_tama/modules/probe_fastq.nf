process PROBE_FASTQ {
    tag "${meta.id}"
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('header_probe.json'), emit: probe
    path 'versions.yml', emit: versions

    script:
    """
    read_input_classification.py probe "${reads}" header_probe.json
    test -s header_probe.json || { echo "Header probe produced no result for ${meta.id}" >&2; exit 1; }
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    printf '{"header_representation":"UNKNOWN","records_sampled":1,"malformed_count":0}\n' > header_probe.json
    printf '"%s":\n    python: 3.11.0\n' '${task.process}' > versions.yml
    """
}
