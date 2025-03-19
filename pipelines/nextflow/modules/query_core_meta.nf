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

process QUERY_CORE_META {
    tag "${core}"
    label 'genomio'
    storeDir "${params.cacheDir}/${core}/meta_data/"

    input:
        val(core)
        path(meta_keys)

    output:
        path("coredb_meta.json"), emit: metadata_json
    
    script:
        """
        genome_metadata_dump --host ${params.host} --port ${params.port} --user ${params.user_r} --database ${core} --metafilter ${meta_keys} --append_db > coredb_meta.json
        """
}