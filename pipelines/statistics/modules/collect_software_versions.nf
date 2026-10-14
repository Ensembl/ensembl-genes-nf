#!/usr/bin/env nextflow
/*
See the NOTICE file distributed with this work for additional information
regarding copyright ownership.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
/*Collect software versions from all processes and merge them into a single file
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