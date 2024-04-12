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

process OMAMER_HOG {

    label 'omamer'
    storeDir "${params.cacheDir}/$gca/omamer/"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    input:
    file fasta
    val gca
    

    output:
    path "proteins.omamer", emit: omamer_file
    val outdir, emit:species_outdir

    script:
    """
    omamer search --db ${params.omamer_database} --query ${fasta} --score sensitive --out proteins.omamer
    
    """
}