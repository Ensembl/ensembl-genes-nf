process INSPECT_BAM_WORKLOAD {
    tag "${meta.id}:bam-workload"
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path(bam), path(bai), path('contig_workload.tsv'), emit: workload
    path 'versions.yml', emit: versions

    script:
    """
    test -r '${bam}' || { echo 'BAM is not readable' >&2; exit 1; }
    test -r '${bai}' || { echo 'BAM index is not readable' >&2; exit 1; }
    samtools quickcheck -v '${bam}'
    samtools idxstats '${bam}' | awk 'BEGIN {OFS="\\t"; print "contig","reference_bases","mapped_reads","eligible","resource_class"} \\
        \$1 != "*" && \$2 > 0 {class=(\$3 >= ${params.shard_contig_reads * 5} ? "very_large" : (\$3 >= ${params.shard_contig_reads} ? "large" : "small")); print \$1,\$2,\$3,"true",class}' > contig_workload.tsv
    test -s contig_workload.tsv || { echo 'BAM has no eligible reference contigs' >&2; exit 1; }
    printf '"%s":\\n    samtools: runtime\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'contig\\treference_bases\\tmapped_reads\\teligible\\tresource_class\\nchrStub\\t4\\t1\\ttrue\\tsmall\\n' > contig_workload.tsv
    printf '"%s":\\n    samtools: stub\\n' '${task.process}' > versions.yml
    """
}
