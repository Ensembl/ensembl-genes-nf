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
    container 'dockerhub.ebi.ac.uk/ensembl_genebuild/ensembl-genes-containers/ensembl-analysis:e52659d38da5af45bb4d8822cef00ea66915bdad'
    storeDir "${params.cacheDir}/${meta.gca}/fasta"
    afterScript "sleep ${params.files_latency}"
    // Needed because of file system latency
    maxForks 20

    input:
    val meta

    output:
    tuple val(meta), path("translations.fa"), emit: fasta_file_output
    path "versions.yml", emit: versions_file

    script:
    def translations_file = "translations.fa"
    """
    if [[ -f "${meta.protein_file}" ]]; then
        echo "Using provided protein file: ${meta.protein_file}"
        cp -L "${meta.protein_file}" ${translations_file}
    else
        dump_translations.pl \
            -host ${params.host} \
            -port ${params.port} \
            -dbname ${meta.dbname} \
            -user ${params.user_r} \
            -file ${translations_file} \
            --species_id ${meta.species_id} \
            ${params.dump_params}
    fi
    # Create versions file - simpler approach
    PERL_VERSION=\$(perl --version | grep -oP 'v\\K[0-9.]+' | head -n1)
    
    echo '"FETCH_PROTEINS":' > versions.yml
    echo "  perl: \$PERL_VERSION" >> versions.yml
    """
}
