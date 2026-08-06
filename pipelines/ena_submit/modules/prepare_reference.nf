process ENA_PREPARE_REFERENCE {
    label 'process_single_long'
    tag { reference_key }
    container 'docker.io/library/python:3.11-slim'

    input:
    tuple val(reference_key), path(reference_fasta), path(assembly_report), path(reference_supplement)

    output:
    tuple val(reference_key), path('reference.fasta'), path('assembly_report.txt'), path('reference.names'), emit: prepared

    script:
    """
    set -euo pipefail
    cat ${reference_fasta} ${reference_supplement} > reference.fasta
    cp ${assembly_report} assembly_report.txt
    awk 'BEGIN { OFS="\\t" } /^>/ { header=substr(\$0, 2); name=header; sub(/[ \\t].*/, "", name); description=header; sub(/^[^ \\t]+[ \\t]*/, "", description); print name, description }' reference.fasta > reference.names
    """
}
