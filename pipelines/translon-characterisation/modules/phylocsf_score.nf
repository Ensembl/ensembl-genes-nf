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
