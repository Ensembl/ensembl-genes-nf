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

include { getMetaValue } from '../utils.nf'

process RUN_STATISTICS {
    label 'fetch_file'
    tag "core_statistics:$gca"
    publishDir "${params.outDir}/$publish_dir/", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    maxForks 20    
    input:
    tuple val(gca), val(dbname)

    output:
    tuple val(gca), val(dbname), path("core_statistics/*.sql")

    script:
    production_name_query = getMetaValue(dbname, "species.production_name")[0]
    production_name = production_name_query ? production_name_query.meta_value.toString() : dbname
    scientific_name_query = getMetaValue(dbname, "species.scientific_name")[0]
    scientific_name = scientific_name_query.meta_value ? scientific_name_query.meta_value.toString().replaceAll("\\s", "_") : dbname
    species=scientific_name.toLowerCase()
    annotation_source_query=getMetaValue(dbname, "species.annotation_source")[0]
    annotation_source = annotation_source_query ? annotation_source_query.meta_value.toString() : "ensembl"
    publish_dir =scientific_name +'/'+gca+'/'+annotation_source
    """
    perl ${params.enscode}/ensembl-genes/src/perl/ensembl/genes/generate_species_homepage_stats.pl \
        -dbname ${dbname} \
        -host ${params.host} \
        -port ${params.port} \
        -production_name ${production_name.trim()} \
        -output_dir core_statistics
    """
    // re write the output to  have a json    
    //gb1-w amphiduros_pacificus_gca949316495v1_core_110_1 <stats_amphiduros_pacificus_gca949316495v1_core_110_1.sql
    

}
