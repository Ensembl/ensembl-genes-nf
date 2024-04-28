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

include { generateMetadataJson } from '../utils.nf'
include { getMetaValue } from '../utils.nf'

process CREATE_STATS_JSON {
    label 'default'
    tag "stats_json:$gca"
    storeDir "${params.outDir}/$publish_dir/core_statistics", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    tuple val(gca), val(core), path(statistics_sql)

    output:
    tuple val(publish_dir), path(statistics_sql)
    path("*.json"), emit: json_files  

    script:
    scientific_name = getMetaValue(core, "species.scientific_name")[0].meta_value.toString().replaceAll("\\s", "_")
    species=scientific_name.toLowerCase()
    publish_dir =scientific_name +'/'+gca+'/'+getMetaValue(core, "species.annotation_source")[0].meta_value.toString()
    gca_string = gca.toLowerCase().replaceAll(/\./, "v").replaceAll(/_/, "")
    json_file_name = [species, gca_string, "core_stats.json"].join("_")
    statistic_file=publish_dir+'/core_statistics/'+statistics_sql.toString()
    //json_file = generateMetadataJson(statistics_sql)
    //output_file= generateMetadataJson(statistic_file, json_file_name)
    // Generate JSON file
    json_files = generateMetadataJson(statistic_file, json_file_name, params.outDir + "/" + publish_dir + "/core_statistics") // Pass outputDir

    
}






