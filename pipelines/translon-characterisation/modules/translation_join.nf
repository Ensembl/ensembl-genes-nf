process TRANSLATION_JOIN {
    label 'process_light'
    container 'oras://community.wave.seqera.io/library/python:3.12-slim'
    tag "${meta.id}"
    publishDir "${params.outdir}/02_translation", mode: 'copy', pattern: '*.reconciled.jsonl'

    input:
    tuple val(meta), path(instances)
    tuple val(verdict_meta), path(verdicts)

    output:
    tuple val(meta), path('*.reconciled.jsonl'), emit: instances
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    """
    translon_characterise.py translation-join --instances ${instances} --verdicts ${verdicts} \\
        --output ${prefix}.reconciled.jsonl ${args}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.reconciled.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
