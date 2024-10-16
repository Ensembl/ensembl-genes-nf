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

process RUN_ENSEMBL_META {

    label 'python'
    conda '../../workflows/bin/python_env.yml'
    tag "$gca"
    publishDir "${params.outDir}/$publish_dir", mode: 'copy'
//    storeDir "${params.cacheDir}/$gca/" 
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
        tuple val(gca), val(dbname)
        output:
        tuple val(gca), val(dbname), path("*.sql")
    
    shell:
        scientific_name = getMetaValue(dbname, "species.scientific_name")[0].meta_value.toString().replaceAll("\\s", "_")
        publish_dir = scientific_name +'/'+gca+'/'+getMetaValue(dbname, "genebuild.annotation_source")[0].meta_value.toString()

        """
        python !{params.enscode}/ensembl-genes/src/python/ensembl/genes/metadata/core_meta_data.py --output_dir !{params.outDir}/$publish_dir/ --db_name !{dbname} --host !{params.host} --port !{params.port}  --team !{params.team}
        ln -s !{params.outDir}/$publish_dir/*.sql 
        """
        //bash mysql -N -u ${params.user} -h ${params.host} -P ${params.port} -D ${dbname} < ${params.cacheDir}/$gca/${dbname}.sql
    

    
}



