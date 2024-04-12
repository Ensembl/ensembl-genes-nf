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

includeConfig '../../../workflows/nextflow.config'
includeConfig '../conf/omark.config'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FETCH_PROTEINS } from '../modules/fetch_proteins.nf'
include { OMAMER_HOG } from '../modules/omark/omamer_hog.nf'
include { OMARK } from '../modules/omark/omark.nf'
include { OMARK_OUTPUT } from '../modules/omark/omark_output.nf'
include { COPY_OUTPUT_TO_ENSEMBL_FTP as COPY_OMARK_OUTPUT } from '../modules/copy_output_to_ensembl_ftp.nf'



include { CLEANING } from '../modules/cleaning.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RUN_OMARK{
    take:                 
    tuple val(dbname),val(db_meta)

    main:
        //
        // MODULE: Get canonical protein from db
        // 
        proteinFile = FETCH_PROTEINS (dbname)
        //
        // MODULE: Get orthologous groups from Omamer db 
        //
        omamerOutput = OMAMER_HOG(proteinFile, db_meta.gca)
        //
        // MODULE: Run Omark
        //        
        omarkOutput = OMARK (omamerOutput)

        omarkSummaryFile = OMARK_OUTPUT(db_meta, omarkOutput, params.project)
        if (params.copyToFtp) {
            COPY_OMARK_OUTPUT(db_meta, omarkSummaryFile)
        }

}


