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
/*DB_METADATA process to fetch core database metadata
Inputs:
- meta: metadata map containing dbname, species_id, taxon_id, gca
Outputs:
- metadata: tuple of meta and metadata.txt file
- versions.yml: versions file containing software versions used
*/
process DB_METADATA {
    label 'python'
    tag "${meta.dbname}"

    input:
    val meta

    output:
    tuple val(meta), path("metadata.txt"), emit: metadata
    path "versions.yml", emit: versions_file

    script:
    """
    if [[ "${meta.taxon_id}" == "UNKNOWN"  &&  "${meta.gca}" == "UNKNOWN" ]]; then
        TAXON_ID=\$(get_meta_value.py \
            --db ${meta.dbname} \
            --key species.taxonomy_id \
            --species-id ${meta.species_id} \
            --host ${params.host} \
            --port ${params.port} \
            --user ${params.user_r})
        GCA=\$(get_meta_value.py \
            --db ${meta.dbname} \
            --key assembly.accession \
            --species-id ${meta.species_id} \
            --host ${params.host} \
            --port ${params.port} \
            --user ${params.user_r}
            )
        PRODUCTION_NAME=\$(get_meta_value.py \
            --db ${meta.dbname} \
            --key species.production_name \
            --species-id ${meta.species_id} \
            --host ${params.host} \
            --port ${params.port} \
            --user ${params.user_r}
            )
    else
        GCA="${meta.gca}";
        TAXON_ID="${meta.taxon_id}";
        PRODUCTION_NAME="${meta.production_name}"
    fi
    # Output metadata to file
    echo "taxon_id=\$TAXON_ID" >> metadata.txt
    echo "gca=\$GCA" >> metadata.txt
    echo "production_name=\$PRODUCTION_NAME" >> metadata.txt
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        get_meta_value.py: \$(get_meta_value.py --version 2>&1 | grep -oP 'version \\K[0-9.]+' || echo "unknown")
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
