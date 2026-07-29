process ENA_REHEADER_BAM {
    label 'process_single_long'
    tag { file_meta.remote_name ?: file.getName() }
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(reference_key), val(meta), val(row), val(file_meta), path(file), path(header)

    output:
    tuple val(reference_key), val(meta), val(row), val(file_meta), path('reheadered.bam'), emit: reheadered

    script:
    """
    set -euo pipefail
    samtools reheader ${header} ${file} > reheadered.bam
    """
}
