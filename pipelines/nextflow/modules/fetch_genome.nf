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


process FETCH_GENOME {
    tag "${organism_name}:${insdc_acc}"
    label 'fetch_file'
    storeDir "${params.cacheDir}/${insdc_acc}/ncbi_dataset/"  // update on protein busco side
    afterScript "sleep ${params.files_latency}"  // Needed because of file system latency
    maxForks 10

    input:
        tuple val(insdc_acc), val(taxonomy_id), val(dbname), 
            val(production_name), val(organism_name), val(annotation_source), val(ortho_db)

    output:
        tuple val(insdc_acc), val(taxonomy_id), val(dbname), 
            val(production_name), val(organism_name), val(annotation_source),
            val(ortho_db), path("*.fna"), emit: genome_fasta
    
    script:
        def outfile = "genome_file.zip"
        """
        curl -X GET "${params.ncbiBaseUrl}/${insdc_acc}/download?include_annotation_type=GENOME_FASTA&hydrated=FULLY_HYDRATED"  -H "Accept: application/zip" --output ${outfile}
        unzip -j $outfile
        """
}
