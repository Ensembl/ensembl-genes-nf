process PREPARE_RIBOCODE_SHARD {
    tag "${meta.id}:${partition.partition_id}"
    label 'process_medium'
    container 'quay.io/biocontainers/ribocode:1.2.15--pyhdc42f0e_1'
    input:
    tuple val(meta), path(bam), path(bai), path(transcriptome_gtf), path(transcriptome_fasta), val(partition)
    output:
    tuple val(meta), path('shard.bam'), path('shard.bam.bai'), path('shard.gtf'), path('shard.fa'), emit: shard
    script:
    def contig = partition.contig
    def start = partition.start as Integer
    def end = partition.end as Integer
    """
    awk -v tx='${contig}' '($0 !~ /^#/ && $0 ~ "transcript_id \\\"" tx "\\\"") || ($0 ~ /^#/)' ${transcriptome_gtf} > shard.gtf
    awk -v tx='${contig}' 'BEGIN{keep=0} /^>/{keep=(substr($0,2)==tx)} keep{print}' ${transcriptome_fasta} > shard.fa
    samtools view -bh ${bam} '${contig}:${start + 1}-${end}' > shard.bam
    samtools index shard.bam
    test -s shard.gtf && test -s shard.fa
    """
    stub:
    """
    touch shard.bam shard.bam.bai
    printf '# stub\n' > shard.gtf
    printf '>stub\nN\n' > shard.fa
    """
}
