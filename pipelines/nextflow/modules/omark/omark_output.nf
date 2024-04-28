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

process OMARK_OUTPUT {
    tag "omark_output:$gca"
    label 'default'
    publishDir "${params.outDir}/$publish_dir/statistics", mode: 'copy'
    //storeDir "${params.outDir}/$publish_dir/statistics/"

    input:
    tuple val(gca), val(dbname), val(publish_dir), path(summary_file), val(omark_dir)

    output:
    tuple val(gca), val(dbname), path("*.txt")

    script:
    println(summary_file)
    scientific_name = getMetaValue(dbname, "species.scientific_name")[0].meta_value.toString().replaceAll("\\s", "_")
    species=scientific_name.toLowerCase()
    gca_string = gca.toLowerCase().replaceAll(/\./, "v").replaceAll(/_/, "")

    summary_name = [species, gca_string, "omark", "proteins_detailed_summary.txt"].join("_")
    summary_file= summary_name

}
