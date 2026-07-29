process ENA_REHEADER_BAM {
    label 'process_single_long'
    tag { file_meta.remote_name ?: file.getName() }
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(meta), val(row), val(file_meta), path(file)
    path reference_fai
    path assembly_report
    path script

    output:
    tuple val(meta), val(row), val(file_meta), path('reheadered.bam'), path('reheadered.bam.bai'), emit: reheadered

    script:
    """
    set -euo pipefail
    python3 ${script} \
        --bam ${file} \
        --reference-fai ${reference_fai} \
        --assembly-report ${assembly_report} \
        --output reheadered.bam
    samtools index -@ ${task.cpus} reheadered.bam reheadered.bam.bai
    samtools quickcheck reheadered.bam
    """
}
