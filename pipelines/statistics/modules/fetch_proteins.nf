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
/*FETCH_PROTEINS process to fetch protein FASTA file from core database
Inputs:
- meta: metadata map containing dbname, species_id, taxon_id, gca
Outputs:
- fasta_file_output: tuple of meta and protein FASTA file
- versions.yml: versions file containing software versions used
*/

process FETCH_PROTEINS {
    tag "${meta.dbname}:protein"
    label 'fetch_file'
    publishDir "${params.cacheDir}/${meta.gca}/fasta", mode: 'copy', pattern: "*.fa"
    publishDir "${params.cacheDir}/${meta.gca}/fasta", mode: 'copy', pattern: "versions.yml"
    afterScript "sleep ${params.files_latency}"
    // Needed because of file system latency
    maxForks 20

    input:
    val meta

    output:
    tuple val(meta), path("translations.fa"), emit: fasta_file_output
    path "versions.yml", emit: versions_file

    script:
    translations_file = "translations.fa"
    """
    if [[ ! -f "${meta.protein_file}" ]]; then
    perl ${params.enscode}/ensembl-analysis/scripts/protein/dump_translations.pl \
        -host ${params.host} \
        -port ${params.port} \
        -dbname ${meta.dbname} \
        -user ${params.user_r} \
        -file ${translations_file} \
        --species_id ${meta.species_id} \
        ${params.dump_params}
    else
    ln -s ${meta.protein_file} ${translations_file}
    fi
    # Create versions file - simpler approach
    PERL_VERSION=\$(perl --version | grep -oP 'v\\K[0-9.]+' | head -n1)
    
    echo '"FETCH_PROTEINS":' > versions.yml
    echo "  perl: \$PERL_VERSION" >> versions.yml
    """
}
