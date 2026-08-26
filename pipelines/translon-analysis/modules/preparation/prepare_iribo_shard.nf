process PREPARE_IRIBO_SHARD {
    tag "${meta.id}:${partition.partition_id}"
    label 'process_high'
    container 'ghcr.io/jackcurragh/translon-iribo:1.0.0'
    input:
    tuple val(meta), path(bam), path(bai), path(gtf), path(fasta), val(partition)
    output:
    tuple val(meta), path('shard.bam'), path('shard.bam.bai'), path('shard.gtf'), path('shard.fa'), emit: shard
    script:
    def contig = partition.contig
    def start = partition.start as Integer
    def end = partition.end as Integer
    """
    samtools view -bh ${bam} '${contig}:${start + 1}-${end}' > shard.bam
    samtools index shard.bam
    # Keep genome coordinates intact: iRibo's GTF uses absolute coordinates.
    cp ${fasta} shard.fa
    awk -v c='${contig}' -v s='${start + 1}' -v e='${end}' '(\$0 !~ /^#/ && \$1 == c && \$5 >= s && \$4 <= e) || (\$0 ~ /^#/)' ${gtf} > shard.gtf
    test -s shard.fa && test -s shard.gtf
    """
    stub:
    """
    touch shard.bam shard.bam.bai
    printf '>stub\nN\n' > shard.fa
    printf '# stub\n' > shard.gtf
    """
}
