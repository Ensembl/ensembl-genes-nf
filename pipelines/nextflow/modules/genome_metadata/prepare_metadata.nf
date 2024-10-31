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

workflow PREPARE_COREDB_METADATA {
    // Generate a meta value from a db metadata file and an ncbi datasets file

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
    label 'default'
    tag "${core}"
    // storeDir "${params.cacheDir}/${core}/meta_data/"

    input:
        val(core)
        path(meta_keys)

    output:
        path("coredb_meta.json"), emit: metadata_json
    

    shell:
        '''
        chmod +x !{projectDir}/bin/meta_data_getter.py
        meta_data_getter.py $(!{params.host} details script) --database_name !{core} --meta_keys_list !{meta_keys}
        '''
}

def meta_from_coredb(json_path) {
    
    //Import meta JSON info from coredb 
    metadata = read_json(json_path)

    // define core_metadata + export for pipeline meta propogation
    return [
        insdc_acc: metadata."assembly.accession",
        taxonomy_id: metadata."species.taxonomy_id",
        dbname: metadata."database_name",
        production_name: metadata."species.production_name",
        organism_name: metadata."species.scientific_name",
        annotation_source: metadata."species.annotation_source",
    ]
}
