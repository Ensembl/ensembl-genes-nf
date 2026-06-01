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

process FETCH_REPEAT_MODEL {
    tag "$meta.gca:repeatmodel"
    label 'fetch_file'
    publishDir "${params.outdir}/${meta.gca}/library/", mode: 'copy'

    input:
    tuple val(meta), path(genome_file)

    output:
    tuple val(meta), path(genome_file), path("${meta.gca}.repeatmodeler.fa"), emit: rep_library_file_output
    path "versions.yml", emit: versions_file


    //tuple val(species_name), val(gca), path(genome_file), path("${gca}.repeatmodeler.fa")
    script:
    """
    # Construct the URL for the repeat model file
    REPEAT_URL="${params.repeats_ftp_base}/${meta.species_name}/${meta.gca}.repeatmodeler.fa"

    # Check if the file exists on the server and download if available
    if curl --silent --fail --output "${meta.gca}.repeatmodeler.fa" "\$REPEAT_URL"; then
        echo "Successfully downloaded RepeatModeler file for ${meta.gca}"
        else
        # Output the GCA if the download is skipped
        echo "Repeat model file not found for ${meta.gca} skipping download"
        echo "No repeatmodeler file available for ${meta.gca}" > "${meta.gca}.repeatmodeler.fa"
    fi

    CURL_VERSION=\$(curl --version 2>&1 | head -n 1 | awk '{print \$2}')
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: \$CURL_VERSION
    END_VERSIONS
    """
}
