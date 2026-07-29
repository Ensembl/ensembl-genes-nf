process ENA_PREPARE_REFERENCE {
    label 'process_single_long'
    tag { reference_key }
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(reference_key), path(reference_fasta), path(assembly_report)

    output:
    tuple val(reference_key), path('reference.fai'), path('assembly_report.txt'), emit: prepared

    script:
    """
    set -euo pipefail
    samtools faidx ${reference_fasta}
    mv ${reference_fasta}.fai reference.fai
    cp ${assembly_report} assembly_report.txt
    """
}
