process ADJUDICATE {
    label 'process_medium'
    container 'oras://community.wave.seqera.io/library/python:3.12-slim'
    tag "${meta.id}"
    publishDir "${params.outdir}/09_adjudication", mode: 'copy', pattern: '*.adjudicated.jsonl'

    input:
    tuple val(meta), path(instances)
    tuple val(cluster_meta), path(clusters)

    output:
    tuple val(meta), path('*.adjudicated.jsonl'), emit: claims
    tuple val(meta), path('*.manifest.json'), emit: manifest
    path 'versions.yml', emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    translon_characterise.py adjudicate --instances ${instances} --clusters ${clusters} --output ${prefix}.adjudicated.jsonl --manifest ${prefix}.manifest.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.adjudicated.jsonl
    printf '{"state":"stub"}\\n' > ${prefix}.manifest.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
