process ENA_INDEX_REFERENCE {
    label 'process_light'
    tag { reference_key }
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(reference_key), path(reference_fasta), path(assembly_report), path(reference_names)

    output:
    tuple val(reference_key), path('reference.fasta'), path('reference.fai'), path('assembly_report.txt'), path('reference.names'), emit: indexed

    script:
    """
    set -euo pipefail
    samtools faidx ${reference_fasta}
    mv ${reference_fasta}.fai reference.fai
    cp ${assembly_report} assembly_report.txt
    cp ${reference_names} reference.names
    """
}
