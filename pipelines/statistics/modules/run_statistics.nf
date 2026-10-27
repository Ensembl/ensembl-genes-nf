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
/*RUN_STATISTICS process to generate core database statistics SQL files
Inputs:
- meta: metadata map containing dbname, species_id, taxon_id, gca
Outputs:
- statistics_output: tuple of meta and generated statistics SQL files
- versions.yml: versions file containing software versions used
*/
process RUN_STATISTICS {
    label 'fetch_file'
    tag "${meta.gca}"
    storeDir "${params.cacheDir}/${meta.gca}/core_statistics/statistics"
    publishDir { "${params.outdir}/${meta.gca}" }, mode: 'copy'
    afterScript "sleep ${params.files_latency}"
    // Needed because of file system latency
    maxForks 20

    input:
    val meta

    output:
    tuple val(meta), path("core_statistics/*.sql"), emit: statistics_output
    path "versions.yml", emit: versions_file

    script:
    """
    perl ${params.enscode}/ensembl-genes/src/perl/ensembl/genes/generate_species_homepage_stats.pl \
        -dbname ${meta.dbname} \
        -host ${params.host} \
        -port ${params.port} \
        -production_name ${meta.production_name} \
        -output_dir core_statistics
    # Create versions file
    PERL_VERSION=\$(perl --version | grep -oP 'v\\K[0-9.]+' | head -n1)
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$PERL_VERSION
    END_VERSIONS
    """
}
