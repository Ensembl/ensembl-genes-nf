process ENA_PREPARE_REFERENCE {
    label 'process_single_long'
    tag { reference_key }
    container 'docker.io/library/python:3.11.13-slim-bookworm'

    input:
    tuple val(reference_key), path(reference_fasta), path(assembly_report), path(reference_supplement)

    output:
    tuple val(reference_key), path('reference.fasta'), path('assembly_report.txt'), path('reference.names'), emit: prepared
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    cat ${reference_fasta} ${reference_supplement} > reference.fasta
    cp ${assembly_report} assembly_report.txt
    awk 'BEGIN { OFS="\\t" } /^>/ { header=substr(\$0, 2); name=header; sub(/[ \\t].*/, "", name); description=header; sub(/^[^ \\t]+[ \\t]*/, "", description); print name, description }' reference.fasta > reference.names
    printf 'ENA_PREPARE_REFERENCE:\n  python: "%s"\n' "\$(python3 --version 2>&1 | awk '{print \$2}')" > versions.yml
    """

    stub:
    """
    touch reference.fasta assembly_report.txt reference.names
    printf 'ENA_PREPARE_REFERENCE:\n  python: "stub"\n' > versions.yml
    """
}
