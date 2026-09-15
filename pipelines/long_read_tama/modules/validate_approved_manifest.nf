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

    script:
    """
    python3 ${validator} validate-approved ${manifest}
    cp ${manifest} validated_approved_run_manifest.tsv
    printf 'status\\tdetail\\nok\\tapproved manifest validated\\n' > approved_manifest_validation.tsv
    """

    stub:
    """
    cp ${manifest} validated_approved_run_manifest.tsv
    printf 'status\\tdetail\\nok\\tstub\\n' > approved_manifest_validation.tsv
    """
}
