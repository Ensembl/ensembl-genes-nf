process indexBAM {

    publishDir output_dir, mode: 'copy'

    input:
        tuple val(sample_id), path(bam)
        path output_dir

    output:
        tuple val(sample_id), path(bam), path("${bam}.bai")

    script:
        """
        samtools index $bam
        """
}
4
process generateBigWig {

    publishDir output_dir, mode: 'copy'

    input:
        tuple val(sample_id), path(bam), path(bai)
        path output_dir

    output:
        tuple val(sample_id), path("${sample_id}_coverage.bw")

    script:
        """
        bamCoverage -b $bam -o ${sample_id}_coverage.bw
        """
}