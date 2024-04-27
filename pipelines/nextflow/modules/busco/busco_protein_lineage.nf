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

process BUSCO_PROTEIN_LINEAGE {
    label 'busco'
    tag "$gca"
    storeDir "${params.cacheDir}/$gca/busco_protein/"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    input:
    val(busco_dataset)
    tuple val(gca), val(dbname), path(translation_file)

    output:
    tuple val(gca), val(dbname), path("protein_output/*.txt")

    script:
    """
    busco -f \
        -i ${translation_file} \
        --mode proteins \
        -l ${busco_dataset.trim()} \
        -c ${task.cpus} \
        --out protein_output \
        --offline \
        --download_path ${params.download_path}
        """
}
