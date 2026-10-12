process ENA_INDEX_BAM {
    label 'process_single_long'
    tag { file_meta.remote_name ?: file.getName() }
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(reference_key), val(meta), val(row), val(file_meta), path(file)

    output:
    tuple val(reference_key), val(meta), val(row), val(file_meta), path('reheadered.bam'), path('reheadered.bam.bai'), emit: indexed
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    samtools index -@ ${task.cpus} ${file} reheadered.bam.bai
    samtools quickcheck ${file}
    printf 'ENA_INDEX_BAM:\n  samtools: "%s"\n' "\$(samtools --version | awk 'NR==1 {print \$2}')" > versions.yml
    """

    stub:
    """
    touch reheadered.bam reheadered.bam.bai
    printf 'ENA_INDEX_BAM:\n  samtools: "stub"\n' > versions.yml
    """
}
