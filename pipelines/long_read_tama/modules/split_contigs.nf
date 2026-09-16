process SPLIT_BAM_BY_CONTIG {
    tag "${meta.id}:contig-shards"
    label 'process_high'
    container 'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1'

    input:
    tuple val(meta), path(bam), path(bai), path(workload)
    path splitter

    output:
    tuple val(meta), path('shards'), emit: shards
    path 'contig_manifest.tsv', emit: manifest
    path 'versions.yml', emit: versions

    script:
    """
    ./${splitter} '${bam}' '${workload}' '${meta.id}' shards contig_manifest.tsv
    printf '"%s":\\n    samtools: runtime\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    mkdir -p shards
    touch shards/${meta.id}.stub.bam shards/${meta.id}.stub.bam.bai
    printf 'contig\\treference_bases\\tmapped_reads\\tresource_class\\tpath\\n' > contig_manifest.tsv
    printf 'chrStub\\t4\\t1\\tsmall\\t${meta.id}.stub.bam\\n' >> contig_manifest.tsv
    printf '"%s":\\n    samtools: stub\\n' '${task.process}' > versions.yml
    """
}
