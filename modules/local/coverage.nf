process indexBAM {

    //publishDir output_dir, mode: 'copy'
    label 'light'

    input:
        tuple val(sample_id), path(bam)
        path output_dir

    output:
        tuple val(sample_id), path(bam), path("${output_dir}/${bam}.bai")

    script:
        """
        samtools index $bam
        """
}
4
process generateBigWig {

    //publishDir output_dir, mode: 'copy'
    label 'light'

    input:
        tuple val(sample_id), path(bam), path(bai)
        path output_dir

    output:
        tuple val(sample_id), path("${output_dir}/${sample_id}_coverage.bw")

    script:
        """
        bamCoverage -b $bam -o ${output_dir}/${sample_id}_coverage.bw
        """
}