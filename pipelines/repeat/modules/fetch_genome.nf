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
    tag "${meta.gca}:genome"
    label 'fetch_file'
    label 'python'
    storeDir "${params.cacheDir}/${meta.gca}/ncbi_dataset"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    maxForks 10

    input:
    val meta

    output:
    tuple val(meta), path("*.fa"), emit: genome_file_output
    path "versions.yml", emit: versions_file

    script:
    """
    if [[ -f "${meta.genome_file}" ]]; then
        echo "Using provided genome file: ${meta.genome_file}"
        cp -L "${meta.genome_file}" genome.fna
    else
        fetch_genome.py \
             --reheader_file \
            --output_dir . \
            --gca ${meta.gca} \
            --ncbi_base ${params.ncbiBaseUrl}

        downloaded_genome=\$(find . -maxdepth 1 -type f -name "*.fa" | head -n 1)

        if [[ -z "\$downloaded_genome" ]]; then
            echo "No genome FASTA found for ${meta.gca}" >&2
            exit 1
        fi
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
