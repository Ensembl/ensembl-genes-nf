process INFLATE_UNIQUE_BAM {
    tag "${meta.merge_group ?: meta.id}:${meta.bam_type}"
    label 'process_long'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    input:
    tuple val(meta), path(bam), path(bai)
    output:
    tuple val(meta), path('inflated.bam'), path('inflated.bam.bai'), emit: inflated
    path 'inflation_manifest.tsv', emit: manifest
    script:
    """
    set -euo pipefail
    samtools view -h ${bam} | awk -v manifest=inflation_manifest.tsv -f ${params.translon_analysis_bin}/inflate_sam.awk | samtools sort -@ ${task.cpus ?: 4} -m 2G -o inflated.bam -
    samtools index -@ ${task.cpus ?: 4} inflated.bam
    """
    stub:
    """
    touch inflated.bam inflated.bam.bai
    printf 'metric\tvalue\nsource_alignments\t0\ninflated_alignments\t0\nsuffixed_alignments\t0\n' > inflation_manifest.tsv
    """
}
