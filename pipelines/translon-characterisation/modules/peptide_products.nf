process PEPTIDE_PRODUCTS {
    label 'process_low'
    container 'oras://community.wave.seqera.io/library/python:3.12-slim'
    tag "${meta.id}"
    input:
    tuple val(meta), path(instances)
    output:
    tuple val(meta), path('*.peptides.jsonl'), emit: peptides
    tuple val(meta), path('*.fanback.jsonl'), emit: fanback
    path 'versions.yml', emit: versions
    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    peptide_products.py --instances ${instances} --peptides ${prefix}.peptides.jsonl --fanback ${prefix}.fanback.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.peptides.jsonl ${prefix}.fanback.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
