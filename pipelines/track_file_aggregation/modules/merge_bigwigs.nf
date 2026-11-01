process MERGE_BIGWIGS {
    tag { "${meta.gca_accession}:${meta.aggregate_id}:${meta.coverage_strand}" }
    label 'process_medium'
    container 'https://depot.galaxyproject.org/singularity/ucsc-bigwigmerge:482--h0b57e2e_0'

    input:
    tuple val(meta), path(bigwigs)

    output:
    tuple val(meta), path("${meta.safe_id}.${meta.coverage_strand}.merged.bedGraph"), emit: bedgraph
    path "${meta.safe_id}.${meta.coverage_strand}.merge.versions.yml", emit: versions

    script:
    """
    printf '%s\\n' *.bw > bigwig_inputs.txt
    bigWigMerge -inList bigwig_inputs.txt ${meta.safe_id}.${meta.coverage_strand}.merged.bedGraph
    test -s ${meta.safe_id}.${meta.coverage_strand}.merged.bedGraph || { echo 'BigWig merge produced no intervals' >&2; exit 1; }
    printf '"%s":\\n    bigWigMerge: 482\\n' '${task.process}' > ${meta.safe_id}.${meta.coverage_strand}.merge.versions.yml
    """

    stub:
    """
    printf 'chr1\\t0\\t100\\t%s\\n' '${meta.source_count}' > ${meta.safe_id}.${meta.coverage_strand}.merged.bedGraph
    printf '"%s":\\n    bigWigMerge: stub\\n' '${task.process}' > ${meta.safe_id}.${meta.coverage_strand}.merge.versions.yml
    """
}
