process PHYLOCSF_SCORE {
    label 'process_medium'
    container params.phylocsf_container
    tag "${meta.id}"

    input:
    tuple val(meta), path(alignments), path(identities)

    output:
    tuple val(meta), path(identities), path('phylocsf.raw.txt'), emit: scored
    path 'versions.yml', emit: versions

    script:
    """
    find ${alignments} -type f -name '*.msa.fasta' -print | sort > alignments.list
    test -s alignments.list
    PhyloCSF ${params.phylocsf_model} --strategy=${params.phylocsf_strategy} --removeRefGaps --files alignments.list > phylocsf.raw.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phylocsf: \$(PhyloCSF --help 2>&1 | head -n1)
        model: ${params.phylocsf_model}
        strategy: ${params.phylocsf_strategy}
    END_VERSIONS
    """

    stub:
    """
    printf 'msas/stub.msa.fasta\t0.0\n' > phylocsf.raw.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phylocsf: stub
        model: ${params.phylocsf_model}
        strategy: ${params.phylocsf_strategy}
    END_VERSIONS
    """
}

process PHYLOCSF_PARSE {
    label 'process_low'
    container 'oras://community.wave.seqera.io/library/python:3.12-slim'
    tag "${meta.id}"
    publishDir "${params.outdir}/06_constraint/phylocsf", mode: 'copy', pattern: '*.phylocsf.axis.jsonl'

    input:
    tuple val(meta), path(identities), path(raw), path(matched_null)

    output:
    tuple val(meta), path('*.phylocsf.axis.jsonl'), emit: axis
    path 'versions.yml', emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    phylocsf_axis.py --identities ${identities} --raw ${raw} --matched-null ${matched_null} --output ${prefix}.phylocsf.axis.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        phylocsf_axis: 1
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.phylocsf.axis.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
        phylocsf_axis: 1
    END_VERSIONS
    """
}
