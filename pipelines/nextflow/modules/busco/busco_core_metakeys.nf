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
    conda '../../workflows/bin/python_env.yml'
    tag "$gca"
    publishDir "${params.outDir}/$publish_dir/", mode: 'copy'
    //storeDir "${params.cacheDir}/$gca/" 
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
        tuple val(gca), val(dbname),path(summary_file)
    
    shell:
        scientific_name = getMetaValue(dbname, "species.scientific_name")[0].meta_value.toString().replaceAll("\\s", "_")
        publish_dir =scientific_name +'/'+gca+'/'+getMetaValue(dbname, "species.annotation_source")[0].meta_value.toString()

        '''
        chmod +x !{projectDir}/bin/busco_metakeys_patch.py
        busco_metakeys_patch.py -db !{dbname} -file !{summary_file} -output_dir "!{params.outDir}/!{publish_dir}/" -host !{params.host} -port !{params.port} -user !{params.user_w}  -password !{params.password} -run_query true
        '''
        // bash mysql -N -u ${params.user} -h ${params.host} -P ${params.port} -D ${dbname} < ${params.cacheDir}/$gca/${dbname}.sql
    
    

    
}



