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

process OMARK {
    label 'omamer'
    tag "$gca"
    
    publishDir "${params.outDir}/$publish_dir/", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    maxForks 15

    input:
    tuple val(gca), val(db), path(omamer_file)
    
    output:
    tuple val(gca), val(db), val(publish_dir), path("omark_output/*.txt"), path("omark_output/*")

    script:
    scientific_name = getMetaValue(db, "species.scientific_name")[0].meta_value.toString().replaceAll("\\s", "_")
    species=scientific_name.toLowerCase()
    publish_dir =scientific_name +'/'+gca+'/'+getMetaValue(db, "species.annotation_source")[0].meta_value.toString()
    
    """
    omark -f ${omamer_file} -d ${params.omamer_database} -o omark_output
    """
}
