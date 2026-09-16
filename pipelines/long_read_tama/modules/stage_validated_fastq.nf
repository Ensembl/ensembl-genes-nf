process STAGE_VALIDATED_FASTQ {
    tag "${meta.id}:${meta.classification}"
    label 'process_light'

    input:
    tuple val(meta), path(reads), path(validation)

    output:
    tuple val(meta), path('reads.fastq.gz'), emit: reads
    path 'versions.yml', emit: versions

    script:
    """
    cp -p "${reads}" reads.fastq.gz
    printf '"%s":\n    staging: complete\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf '@stub/ccs\nACGT\n+\n!!!!\n' | gzip -c > reads.fastq.gz
    printf '"%s":\n    staging: complete\n' '${task.process}' > versions.yml
    """
}
