process ENA_BUILD_BAM_HEADER {
    label 'process_light'
    tag { file_meta.remote_name ?: file.getName() }
    container 'docker.io/library/python:3.11.13-slim-bookworm'

    input:
    tuple val(reference_key), val(meta), val(row), val(file_meta), path(file), path(header), path(reference_fai), path(assembly_report), path(reference_names), path(script)

    output:
    tuple val(reference_key), val(meta), val(row), val(file_meta), path('source.bam'), path('reheader.header'), emit: built
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    python3 ${script} \
        --header ${header} \
        --reference-fai ${reference_fai} \
        --assembly-report ${assembly_report} \
        --reference-names ${reference_names} \
        --output-header reheader.header
    printf 'ENA_BUILD_BAM_HEADER:\n  python: "%s"\n' "\$(python3 --version 2>&1 | awk '{print \$2}')" > versions.yml
    """

    stub:
    """
    touch source.bam reheader.header
    printf 'ENA_BUILD_BAM_HEADER:\n  python: "stub"\n' > versions.yml
    """
}
