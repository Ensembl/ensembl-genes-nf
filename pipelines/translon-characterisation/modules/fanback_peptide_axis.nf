process FANBACK_PEPTIDE_AXIS {
    label 'process_light'
    container 'oras://community.wave.seqera.io/library/python:3.12-slim'
    tag "${meta.id}"
    input:
    tuple val(meta), path(instances)
    tuple val(peptide_meta), path(peptide_axes)
    output:
    tuple val(meta), path('*.peptide_uniqueness.axis.jsonl'), emit: axis
    path 'versions.yml', emit: versions
    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    fanback_peptide_axis.py --instances ${instances} --peptides ${peptide_axes} --output ${prefix}.peptide_uniqueness.axis.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.peptide_uniqueness.axis.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
