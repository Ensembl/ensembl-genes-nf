process ENA_BUILD_BAM_HEADER {
    label 'process_light'
    tag { file_meta.remote_name ?: file.getName() }
    container 'docker.io/library/python:3.11-slim'

    input:
    tuple val(reference_key), val(meta), val(row), val(file_meta), path(file), path(header), path(reference_fai), path(assembly_report), path(reference_names), path(script)

    output:
    tuple val(reference_key), val(meta), val(row), val(file_meta), path('source.bam'), path('reheader.header'), emit: built

    script:
    """
    set -euo pipefail
    python3 ${script} \
        --header ${header} \
        --reference-fai ${reference_fai} \
        --assembly-report ${assembly_report} \
        --reference-names ${reference_names} \
        --output-header reheader.header
    """
}
