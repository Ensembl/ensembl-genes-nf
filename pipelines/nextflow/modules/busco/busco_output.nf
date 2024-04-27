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

process BUSCO_OUTPUT {
    label 'default'
    tag "busco_output:$gca"
    publishDir "${params.outDir}/$publish_dir/statistics", mode: 'copy'
    
    input:
    val(datatype)
    tuple val(core), val(gca), path(summary_file) 


    output:
    val publish_dir, emit: output_dir
    path("*_short_summary.txt"), emit:renamed_summary_file
    
    script:
    if (core=='core'){
         publish_dir =gca 
         species="species"
    }else{
         scientific_name = getMetaValue(core, "species.scientific_name")[0].meta_value.toString().replaceAll("\\s", "_")
         species=scientific_name.toLowerCase()
         publish_dir =scientific_name +'/'+gca+'/'+getMetaValue(core, "species.annotation_source")[0].meta_value.toString()
    }
    
    def name = ""
    if (datatype == "genome") {
        name = "genome_busco"
    } else if (datatype == "protein") {
            name = "protein_busco"
    }
    if (params.project == 'brc') {
        summary_name = [name, "short_summary.txt"].join("_")
    }
    else {
        gca_string = gca.toLowerCase().replaceAll(/\./, "v").replaceAll(/_/, "")
        summary_name = [species, gca_string, name, "short_summary.txt"].join("_")
    } 
    """
    sed '/Summarized benchmarking in BUSCO notation for file/d' $summary_file > $summary_name
    """
}
