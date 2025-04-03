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
    tag "${organism_name}:${insdc_acc}"
    storeDir "${params.cacheDir}/${insdc_acc}/omamer/"
    container "${params.omark_singularity_path}"
    
    input:
        tuple val(insdc_acc), val(taxonomy_id), val(dbname), 
            val(production_name), val(organism_name), val(annotation_source),
            val(ortho_db), path (translation_seqs)

    output:
        tuple val(insdc_acc), val(taxonomy_id), val(dbname), 
            val(production_name), val(organism_name), val(annotation_source), path("proteins.omamer")

    script:
        """
        omamer search --db ${params.omamer_database} --query ${translation_seqs} --score sensitive --out proteins.omamer
        """
}
