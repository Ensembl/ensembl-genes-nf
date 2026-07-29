process ENA_INDEX_CRAM {
    label 'process_single_long'
    tag { file_meta.remote_name ?: file.getName() }
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(meta), val(row), val(file_meta), path(file)

    output:
    tuple val(meta), val(row), val(file_meta), path('converted.cram'), path('converted.cram.crai'), emit: indexed

    script:
    """
    set -euo pipefail
    samtools index -@ ${task.cpus} ${file} converted.cram.crai
    samtools quickcheck ${file}
    """
}
