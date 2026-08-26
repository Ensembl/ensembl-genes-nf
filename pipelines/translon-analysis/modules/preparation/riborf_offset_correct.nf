process PREPARE_RIBORF_READS {
    tag "${meta.id}:${meta.shard_id ?: 'all'}"
    label 'process_low'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'ghcr.io/jackcurragh/translon-riborf:1.0.0'
    input:
    tuple val(meta), path(bam), path(bai), path(genepred), path(bed12), path(offsets)
    output:
    tuple val(meta), path('reads.sam'), path(genepred), path(bed12), emit: riborf
    script:
    """
    samtools view -h ${bam} > uncorrected.sam
    perl \$RIBORF_HOME/offsetCorrect.pl -r uncorrected.sam -p ${offsets} -o reads.sam
    printf '"%s":\n    samtools: \$(samtools --version | head -n1)\n' '${task.process}' > versions.yml
    """
    stub:
    """
    touch reads.sam
    printf '"stub":\n    samtools: stub\n' > versions.yml
    """
}
