// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Utilities
include { read_json } from '../utils.nf'

workflow PREPARE_METADATA {
    // Generate a meta object from a database meta table and dump in JSON format

    take:
        input_cores
        core_metakeys // path -> input meta keys TXT

    emit:
        core_metadata
    
    main:
        _QUERY_CORE_META(input_cores, core_metakeys)
            .map{ meta_json_file -> meta_from_coredb(meta_json_file) }
            .set { core_metadata }
}


process _QUERY_CORE_META {
    tag "${core}"
    label 'genomio'
    storeDir "${params.cacheDir}/${core}/meta_data/"
    afterScript "sleep ${params.files_latency}"

    input:
        val(core)
        path(meta_keys)

    output:
        path("coredb_meta.json"), emit: metadata_json
    
    shell:
        '''
        genome_metadata_dump --host !{params.host} --port !{params.port} --user !{params.user_r} --database !{core} --metafilter !{meta_keys} --append_db > coredb_meta.json
        '''
}

def meta_from_coredb(json_path) {
    
    //Import meta JSON info from coredb 
    metadata = read_json(json_path)

    // define core_metadata + export for pipeline meta propogation
    return [
        insdc_acc: metadata.get("assembly").get("accession"),
        taxonomy_id: metadata.get("species").get("taxonomy_id"),
        dbname: metadata.get("database").get("name"),
        production_name: metadata.get("species").get("production_name"),
        organism_name: metadata.get("species").get("scientific_name"),
        annotation_source: metadata.get("species").get("annotation_source"),
    ]
}
