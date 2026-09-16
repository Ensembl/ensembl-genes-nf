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
    inspect_bam_workload.sh '${bam}' '${bai}' contig_workload.tsv ${params.shard_contig_reads}
    printf '"%s":\\n    samtools: runtime\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'contig\\treference_bases\\tmapped_reads\\teligible\\tresource_class\\nchrStub\\t4\\t1\\ttrue\\tsmall\\n' > contig_workload.tsv
    printf '"%s":\\n    samtools: stub\\n' '${task.process}' > versions.yml
    """
}
