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

process OMARK {
    label 'omamer'
    tag "$meta.gca"

    publishDir "${params.outdir}/$meta.gca", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    maxForks 15

    input:
    tuple val(meta), path(omamer_file)

    output:
    tuple val(meta), path("omark_output/*.txt"), path("omark_output/*"), emit: omark_output
    path "versions.yml", emit: versions_file

    script:
    """
    omark -f ${omamer_file} -d ${params.omamer_database} -o omark_output
    # Create versions file
    # OMARK_VERSION=\$(omark --version 2>&1 | grep -oP 'OMark\\s+v?\\K[0-9.]+' || echo "unknown")
    OMARK_VERSION=\$(pip show omark | grep -i '^Version:' | awk '{print \$2}')    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omark: \$OMARK_VERSION
    END_VERSIONS
    """
}
