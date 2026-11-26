#!/usr/bin/env nextflow
/*
Collect software versions from all processes and merge them into a single file
Inputs:
- versions_*.yml: versions files from all processes
Outputs:- software_versions.yml: merged versions file containing software versions used in the pipeline

*/
// Single process to merge all versions
process COLLECT_SOFTWARE_VERSIONS {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path 'versions_*.yml'

    output:
    path "software_versions.yml"

    script:
    """
    cat versions_*.yml > software_versions.yml
    """
}