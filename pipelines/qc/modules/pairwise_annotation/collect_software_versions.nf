process COLLECT_SOFTWARE_VERSIONS {
    tag 'pairwise_qc_versions'
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', pattern: 'software_versions.yml'

    input:
    path version_files, stageAs: 'versions_*.yml'

    output:
    path 'software_versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    cat versions_*.yml > software_versions.yml
    """

    stub:
    """
    printf '# process versions\\n' > software_versions.yml
    """
}
