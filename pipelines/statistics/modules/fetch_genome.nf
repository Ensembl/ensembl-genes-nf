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

/*FETCH_GENOME process to fetch genome FASTA file from NCBI using GCA accession
Inputs:
- meta: metadata map containing gca
Outputs:
- fasta_file_output: tuple of meta and genome FASTA file
- versions.yml: versions file containing software versions used
*/
process FETCH_GENOME {
    tag "${meta.gca}:genome"
    label 'fetch_file'
    label 'python'


    input:
    val meta

    output:
    tuple val(meta), path("genome.fna"), emit: fasta_file_output
    path "versions.yml", emit: versions_file

    script:
    def cache_dir = "${params.cacheDir}/${meta.gca}/ncbi_dataset"
    """
    mkdir -p "${cache_dir}"

    if [[ -f "${meta.genome_file}" ]]; then
        ln -s "${meta.genome_file}" genome.fna
    else
        cached_genome=\$(find "${cache_dir}" -maxdepth 1 -type f -name "*.fna" | head -n 1)

        if [[ -z "\$cached_genome" ]]; then
            fetch_genome.py \
                --output_dir "${cache_dir}" \
                --gca ${meta.gca} \
                --ncbi_base ${params.ncbiBaseUrl}

            cached_genome=\$(find "${cache_dir}" -maxdepth 1 -type f -name "*.fna" | head -n 1)
        fi

        if [[ -z "\$cached_genome" ]]; then
            echo "No genome FASTA found for ${meta.gca}" >&2
            exit 1
        fi

        ln -s "\$cached_genome" genome.fna
    fi
    
    # Create versions file
    PYTHON_VERSION=\$(python --version 2>&1 | awk '{print \$2}')
    FETCH_GENOME_VERSION=\$(fetch_genome.py  --version 2>&1 || echo "unknown")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fetch_genome.py: \$FETCH_GENOME_VERSION
        python: \$PYTHON_VERSION
    END_VERSIONS
    
    """
}
