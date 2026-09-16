process VALIDATE_FASTQ {
    tag "${meta.id}:${meta.classification}"
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('full_validation.json'), emit: validation
    path 'versions.yml', emit: versions

    script:
    def expected = meta.expected_header_representation ?: 'UNKNOWN'
    """
    read_input_classification.py validate-fastq "${reads}" "${expected}" full_validation.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    printf '{"records":1,"representation":"UNKNOWN","distinct_ids":1,"distinct_molecules":1}\n' > full_validation.json
    printf '"%s":\n    python: 3.11.0\n' '${task.process}' > versions.yml
    """
}
