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

nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COPY_OUTPUT_TO_ENSEMBL_FTP as COPY_OMARK_OUTPUT } from '../modules/copy_output_to_ensembl_ftp.nf'
include { FETCH_PROTEINS } from '../modules/fetch_proteins.nf'
include { OMARK } from '../modules/omark/omark.nf'
include { OMAMER_HOG } from '../modules/omark/omamer_hog.nf'
include { OMARK_OUTPUT } from '../modules/omark/omark_output.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RUN_OMARK{
    take:                 
        db_meta

    main:
        // Need to alter hashmap to add dummy variable to fetch proteins module
        db_meta
            .map{ it -> [
                insdc_acc: it['insdc_acc'], taxonomy_id: it['taxonomy_id'], dbname: it['dbname'], production_name: it['production_name'],
                organism_name: it['organism_name'], annotation_source: it['annotation_source'], ortho_db: "NoOrthoDB"
            ]}
            .set{ amended_db_meta }

        // MODULE: Get canonical protein from db
        def proteinData = FETCH_PROTEINS(amended_db_meta)

        // MODULE: Get orthologous groups from Omamer db 
        def omamerData = OMAMER_HOG(proteinData)

        // MODULE: Run Omark
        def omarkOutput = OMARK(omamerData)

        // Output Omark results
        def omarkSummaryOutput = OMARK_OUTPUT(omarkOutput)


}


