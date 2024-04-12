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
    label 'fetch_input'
    tag "$db_meta.species:$db_meta.gca"
    storeDir "${params.cacheDir}/$gca/omark_output"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    input:
    val db_meta

    output:
    path("proteins_detailed_summary.txt"), emit: summary_file

    script:
    """
    perl ${params.enscode}/core_meta_updates/scripts/stats/generate_species_homepage_stats.pl \
        -dbname ${db_meta.name} \
        -host ${params.host} \
        -port ${params.port} \
        -production_name ${db_meta.production_name}
    // re write the output to  have a json    
    gb1-w amphiduros_pacificus_gca949316495v1_core_110_1 <stats_amphiduros_pacificus_gca949316495v1_core_110_1.sql
    
    """
}