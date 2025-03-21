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

process BUSCO_DATASET {
    label 'python'
    conda "${projectDir}/bin/environment.yml"
    tag "${organism_name} [tax_id:${taxonomy_id}]"
    // storeDir "${params.cacheDir}/${dbname}/meta_data/"

    input:
        tuple val(insdc_acc), val(taxonomy_id), val(dbname),
            val(production_name), val(organism_name), val(annotation_source)

    output:
        tuple val(insdc_acc), val(taxonomy_id), val(dbname),
            val(production_name), val(organism_name), val(annotation_source),
            env(orthodb), emit: clade_dataset
        path("orthodb_set.txt")
    
    script:
        output = 'orthodb_set.txt'
        orthodb = ''
        """
        clade_selector.py -d $params.busco_datasets_file -t $taxonomy_id > $output
        orthodb=`cat ${output}`
        """
}
