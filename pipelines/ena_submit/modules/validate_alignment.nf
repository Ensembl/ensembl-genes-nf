process ENA_VALIDATE_ALIGNMENT {
    label 'process_single_long'
    tag { file_meta.remote_name ?: file.getName() }
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(reference_key), val(meta), val(row), val(file_meta), path(file), path(reference_fai), path(assembly_report), path(reference_names), path(script)

    output:
    tuple val(meta), val(row), val(file_meta), path('validated.alignment'), emit: validated
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    samtools quickcheck ${file}
    samtools view -H ${file} > input.header
    python3 ${script} \\
        --header input.header \\
        --reference-fai ${reference_fai} \\
        --assembly-report ${assembly_report} \\
        --reference-names ${reference_names} \\
        --output-header validated.header
    cp ${file} validated.alignment
    printf 'ENA_VALIDATE_ALIGNMENT:\\n  samtools: "%s"\\n' "\\$(samtools --version | awk 'NR==1 {print \\$2}')" > versions.yml
    """

    stub:
    """
    touch validated.alignment
    printf 'ENA_VALIDATE_ALIGNMENT:\\n  samtools: "stub"\\n' > versions.yml
    """
}
