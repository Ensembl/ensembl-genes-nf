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
/*This process fetches the RepeatModeler library file for a given GCA accession. 
If the file is available at the specified URL, it downloads the file 
and saves it with the name "<GCA>.repeatmodeler.fa". 
If the file is not available, it outputs a message indicating that 
the download was skipped and creates an empty file with the same name. 
The fetched library file is saved in the "library" directory under 
the output directory for the given GCA accession.*/

process FETCH_REPEAT_MODEL {
    tag "$meta.gca:repeatmodel"
    label 'fetch_file'
    publishDir "${params.outdir}/${meta.gca}/library/", mode: 'copy'

    input:
    val(meta)

    output:
    tuple val(meta), path("${meta.gca}.repeatmodeler.fa"), emit: rep_library_file_output
    path "versions.yml", emit: versions_file

    script:
    """
    # Construct the URL for the repeat model file
    REPEAT_URL="${params.repeats_ftp_base}/${meta.species_name}/${meta.gca}.repeatmodeler.fa"
    echo "\$REPEAT_URL"
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
