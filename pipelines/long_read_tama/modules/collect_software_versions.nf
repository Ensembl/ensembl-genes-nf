process COLLECT_LONG_READ_SOFTWARE_VERSIONS {
    tag 'software-versions'
    label 'process_light'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path 'versions_*.yml'

    output:
    path 'software_versions.yml', emit: versions

    script:
    """
    cat versions_*.yml > software_versions.yml
    """

    stub:
    """
    cat versions_*.yml > software_versions.yml
    """
}
