process ENA_INDEX_BAM {
    label 'process_single_long'
    tag { file_meta.remote_name ?: file.getName() }
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(reference_key), val(meta), val(row), val(file_meta), path(file)

    output:
    tuple val(reference_key), val(meta), val(row), val(file_meta), path('reheadered.bam'), path('reheadered.bam.bai'), emit: indexed

    script:
    """
    set -euo pipefail
    samtools index -@ ${task.cpus} ${file} reheadered.bam.bai
    samtools quickcheck ${file}
    """
}
