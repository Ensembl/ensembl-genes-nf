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

process OMAMER_HOG {
    maxForks 15
    label 'omamer'
    tag "${meta.gca}"
    storeDir "${params.cacheDir}/${meta.gca}/omamer/"
    afterScript "sleep ${params.files_latency}"

    input:
    tuple val(meta), path(translation_file)

    output:
    tuple val(meta), path("proteins.omamer"), emit: omamer_hog_output
    path "versions.yml", emit: versions_file

    script:
    """
        omamer search --db ${params.omamer_database} --query ${translation_file}  --out proteins.omamer 

        # Create versions file
        #OMAMER_VERSION=\$(omamer --version 2>&1 | grep -oP 'OMAmer\\s+v?\\K[0-9.]+' || echo "unknown")
    
        OMAMER_VERSION=\$(pip show omamer | grep -i '^Version:' | awk '{print \$2}')

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            omamer: \$OMAMER_VERSION
        END_VERSIONS
        """
}
