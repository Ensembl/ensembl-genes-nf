process TRIM_CAGE_READS {
    tag "${meta.sample_id}:${meta.protocol}"
    label 'cage_trim'
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.trim.fastq*'), emit: reads

    script:
    def input_args = reads.size() == 2 ? "-o ${reads[0].simpleName}.trim.fastq -p ${reads[1].simpleName}.trim.fastq ${reads[0]} ${reads[1]}" : "-o ${reads[0].simpleName}.trim.fastq ${reads[0]}"
    """
    cutadapt ${input_args} ${meta.trim_args ?: ''} > ${meta.sample_id}.cutadapt.log
    """

    stub:
    """
    for read in ${reads.join(' ')}; do cp "\$read" "\${read%.fastq}.trim.fastq"; done
    """
}
