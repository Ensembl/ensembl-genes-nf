process BEDGRAPH_TO_BIGWIG {
    tag { "${meta.gca_accession}:${meta.aggregate_id}:${meta.coverage_strand}" }
    label 'process_medium'
    container 'https://depot.galaxyproject.org/singularity/ucsc-bedgraphtobigwig:482--hdc0a859_0'

    input:
    tuple val(meta), path(bedgraph), path(chrom_sizes)

    output:
    tuple val(meta), path("${meta.safe_id}.${meta.coverage_strand}.bw"), emit: bigwig
    path "${meta.safe_id}.${meta.coverage_strand}.bigwig.versions.yml", emit: versions

    script:
    """
    bedGraphToBigWig ${bedgraph} ${chrom_sizes} ${meta.safe_id}.${meta.coverage_strand}.bw
    test -s ${meta.safe_id}.${meta.coverage_strand}.bw || { echo 'Aggregated BigWig is empty' >&2; exit 1; }
    printf '"%s":\\n    bedGraphToBigWig: 482\\n' '${task.process}' > ${meta.safe_id}.${meta.coverage_strand}.bigwig.versions.yml
    """

    stub:
    """
    printf 'stub\\n' > ${meta.safe_id}.${meta.coverage_strand}.bw
    printf '"%s":\\n    bedGraphToBigWig: stub\\n' '${task.process}' > ${meta.safe_id}.${meta.coverage_strand}.bigwig.versions.yml
    """
}
