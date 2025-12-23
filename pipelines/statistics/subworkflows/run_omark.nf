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

nextflow.enable.dsl=2


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FETCH_PROTEINS } from '../modules/fetch_proteins.nf'
include { OMAMER_HOG } from '../modules/omamer_hog.nf'
include { OMARK } from '../modules/omark.nf'
//include { OMARK_OUTPUT } from '../modules/omark/omark_output.nf'
//include { COPY_OUTPUT_TO_ENSEMBL_FTP as COPY_OMARK_OUTPUT } from '../modules/copy_output_to_ensembl_ftp.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RUN_OMARK{
    take:                 
    csvFile

    main:
    // Validate required parameters
    if (params.fetch == null) {
        error "params.fetch must be defined as true or false"
    }
    // Read data from the CSV file, split it, and map each row to extract GCA and taxon values
        data = Channel.fromPath(csvFile, type: 'file', checkIfExists: true)
                .splitCsv(sep:',', header:true)
                .map { row -> 
                    [gca:row.get('gca'), 
                    core:row.get('core'),
                    species_id:row.get('species_id')
                    ]
                    }
        def db_meta=data
        // MODULE: Get canonical protein from db
        // 
        def proteinData = FETCH_PROTEINS (db_meta).output.protein_file_output
        //
        // MODULE: Get orthologous groups from Omamer db 
        //
        def omamerData = OMAMER_HOG(proteinData).output.omamer_hog_output
        //
        // MODULE: Run Omark
        //        
        OMARK (omamerData).output.omark_output

        //def (omarkSummaryOutput) = OMARK_OUTPUT(omarkOutput)
        //if (params.copyToFtp) {
        //    COPY_OMARK_OUTPUT(omarkSummaryOutput)
        //}

}


