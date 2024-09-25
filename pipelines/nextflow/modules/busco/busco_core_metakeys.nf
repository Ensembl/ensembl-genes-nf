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

process BUSCO_CORE_METAKEYS {

    label 'python'
    tag "$gca"
    storeDir "${params.cacheDir}/$gca/" 
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    tuple val(gca), val(dbname),path(summary_file)
    
    script:
    """
    busco_metakeys_patch.py -db ${dbname} -file ${summary_file}
    """

    
}



