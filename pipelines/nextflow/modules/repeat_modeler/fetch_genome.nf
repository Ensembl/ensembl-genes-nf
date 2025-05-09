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

process FETCH_GENOME {
    tag "$gca:genome"
    label 'fetch_file'
    publishDir "${params.outDir}/${gca}", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    maxForks 10

    input:
    tuple val(species), val(gca)

    output:
    tuple val(species),val(gca),path("*.fasta")

    script:
    """
    set -e
    curl -X GET "${params.ncbiBaseUrl}/${gca}/download?include_annotation_type=GENOME_FASTA&hydrated=FULLY_HYDRATED" -H "Accept: application/zip" --output genome_file.zip
    unzip -j genome_file.zip
    for file in *.fna; do
        mv "\$file" "\${file%.fna}.fasta"
    done
    """
}
