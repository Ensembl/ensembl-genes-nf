process COLLECT_SOFTWARE_VERSIONS {
    label 'process_light'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path 'versions_*.yml'

    output:
    path 'software_versions.yml'

    script:
    """
    set -euo pipefail
    cat versions_*.yml > software_versions.yml
    """

    stub:
    """
    cat versions_*.yml > software_versions.yml
    """
}
