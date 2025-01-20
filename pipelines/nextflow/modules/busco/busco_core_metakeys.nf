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

process BUSCO_CORE_METAKEYS {

    label 'python'
    tag "$gca"
    publishDir "${params.outDir}/$publish_dir/", mode: 'copy'
    //storeDir "${params.cacheDir}/$gca/" 
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    tuple val(gca), val(dbname),path(summary_file)
    
    script:
    scientific_name_query = getMetaValue(dbname, "species.scientific_name")[0]
    scientific_name = scientific_name_query.meta_value ? scientific_name_query.meta_value.toString().replaceAll("\\s", "_") : dbname
    species=scientific_name.toLowerCase()
    annotation_source_query=getMetaValue(dbname, "species.annotation_source")[0]
    annotation_source = annotation_source_query ? annotation_source_query.meta_value.toString() : "ensembl"
    publish_dir =scientific_name +'/'+gca+'/'+annotation_source

    """
    # Check if Python dependencies are installed
    # Read each line in the requirements file
    while read -r package; do \\
    if ! pip show -q "\$package" &>/dev/null; then 
        echo "\$package is not installed" 
        pip install "\$package"
    else
        echo "\$package is already installed"
    fi
    done < ${projectDir}/bin/requirements.txt

    chmod +x $projectDir/bin/busco_metakeys_patch.py
    busco_metakeys_patch.py -db ${dbname} -file ${summary_file} -output_dir "${params.outDir}/$publish_dir/" -host ${params.host} -port ${params.port} -user ${params.user}  -password ${params.password} -run_query true
    """
    //bash mysql -N -u ${params.user} -h ${params.host} -P ${params.port} -D ${dbname} < ${params.cacheDir}/$gca/${dbname}.sql
    

    
}



