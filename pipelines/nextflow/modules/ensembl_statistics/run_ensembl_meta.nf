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

process RUN_ENSEMBL_META {

    label 'python'
    conda '../../workflows/bin/python_env.yml'
    tag "${insdc_acc}"
    publishDir "${params.outDir}/${publish_dir_name}", mode: 'copy'
    // afterScript "sleep ${params.files_latency}"  // Needed because of file system latency
    // storeDir "${params.cacheDir}/$insdc_acc/" 

    input:
        tuple val(insdc_acc), val(taxonomy_id), val(dbname), 
            val(production_name), val(organism_name), val(annotation_source)

    output:
                tuple val(insdc_acc), val(taxonomy_id), val(dbname), 
            val(production_name), val(organism_name), val(annotation_source), path("*.sql")
    
    shell:
        formated_sci_name = organism_name.replaceAll("\\s", "_")
        publish_dir_name = formated_sci_name + '/' + insdc_acc + '/' + annotation_source

        '''
        python !{params.enscode}/ensembl-genes/src/python/ensembl/genes/metadata/core_meta_data.py --output_dir !{params.outDir}/!{publish_dir_name}/ --db_name !{dbname} --host !{params.host} --port !{params.port}  --team !{params.team}
        ln -s !{params.outDir}/!{publish_dir_name}/*.sql 
        '''
        //bash mysql -N -u ${params.user} -h ${params.host} -P ${params.port} -D ${dbname} < ${params.cacheDir}/$gca/${dbname}.sql
    

    
}



