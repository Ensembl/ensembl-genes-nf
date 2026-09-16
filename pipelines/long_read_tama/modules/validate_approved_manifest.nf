process VALIDATE_APPROVED_LONG_READ_MANIFEST {
    tag 'approved-manifest'
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    path manifest
    path validator

    output:
    path 'validated_approved_run_manifest.tsv', emit: manifest
    path 'approved_manifest_validation.tsv', emit: report
    path 'versions.yml', emit: versions

    script:
    """
    python3 ${validator} validate-approved ${manifest}
    test \$(awk 'NR > 1 && NF > 0 {n++} END {print n+0}' ${manifest}) -gt 0 || { echo "Approved manifest is empty: ${manifest}" >&2; exit 1; }
    cp ${manifest} validated_approved_run_manifest.tsv
    printf 'status\\tdetail\\nok\\tapproved manifest validated\\n' > approved_manifest_validation.tsv
    printf '"%s":\\n    python: runtime\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    cp ${manifest} validated_approved_run_manifest.tsv
    printf 'status\\tdetail\\nok\\tstub\\n' > approved_manifest_validation.tsv
    printf '"%s":\\n    python: stub\\n' '${task.process}' > versions.yml
    """
}
