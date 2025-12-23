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

process RUN_STATISTICS {
    label 'fetch_file'
    tag "$meta.gca"
    publishDir "${params.cacheDir}/$meta.gca/core_statistics", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    maxForks 20    
    input:
    val(meta)

    output:
    tuple val(meta), path("core_statistics/*.sql"), emit: statistics_output
    path "versions.yml", emit: versions_file

    script:

    """
    PRODUCTION_NAME=\$(python utils.py \
    --db ${meta.core} \
    --key species.production_name \
    --species-id ${meta.species_id} \
    --host ${params.host} \
    --port ${params.port} \
    --user ${params.user_r}
    )

    perl ${params.enscode}/ensembl-genes/src/perl/ensembl/genes/generate_species_homepage_stats.pl \
        -dbname ${meta.core} \
        -host ${params.host} \
        -port ${params.port} \
        -production_name \$PRODUCTION_NAME \
        -output_dir core_statistics
    """

}
