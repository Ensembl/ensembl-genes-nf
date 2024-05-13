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

//includeConfig '../../../workflows/nextflow.config'
//includeConfig '../conf/busco.config'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BUSCO_DATASET } from '../modules/busco/busco_dataset.nf'
include { FETCH_GENOME } from '../modules/fetch_genome.nf'
include { FETCH_PROTEINS } from '../modules/fetch_proteins.nf'
include { BUSCO_GENOME_LINEAGE } from '../modules/busco/busco_genome_lineage.nf'
include { BUSCO_PROTEIN_LINEAGE } from '../modules/busco/busco_protein_lineage.nf'
include { BUSCO_OUTPUT as BUSCO_GENOME_OUTPUT } from '../modules/busco/busco_output.nf'
include { BUSCO_OUTPUT as BUSCO_PROTEIN_OUTPUT } from '../modules/busco/busco_output.nf'
include { COPY_OUTPUT_TO_ENSEMBL_FTP as COPY_GENOME_OUTPUT } from '../modules/copy_output_to_ensembl_ftp.nf'
include { COPY_OUTPUT_TO_ENSEMBL_FTP as COPY_PROTEIN_OUTPUT } from '../modules/copy_output_to_ensembl_ftp.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow RUN_BUSCO{
    take:                 
    db_meta
    busco_mode
    copyToFtp

    main:
    // Get the closest Busco dataset from the taxonomy classification stored in db meta table
    def db_meta1=db_meta
    db_meta1.flatten().view { d -> "GCA: ${d.gca}, Taxon ID: ${d.taxon_id}, Core name: ${d.core}"}
    
    def (dataset, db) = BUSCO_DATASET(db_meta.flatten()) 
    // Run Busco in genome mode
    if (busco_mode.contains('genome')) {
        def output_typeG = "genome"
        def genomeData = FETCH_GENOME(db)
        def buscoGenomeOutput = BUSCO_GENOME_LINEAGE(dataset, genomeData)
        def buscoGenomeSummaryOutput = BUSCO_GENOME_OUTPUT(output_typeG,buscoGenomeOutput)
        if (params.copyToFtp) {
            COPY_GENOME_OUTPUT(buscoGenomeSummaryOutput)
        }
    }
    
    // Run Busco in protein mode
    if (busco_mode.contains('protein')) {
        def output_typeP = "protein"
        def proteinData = FETCH_PROTEINS (db)
        def buscoProteinOutput = BUSCO_PROTEIN_LINEAGE(dataset,proteinData)
        def (buscoProteinSummaryOutput) = BUSCO_PROTEIN_OUTPUT(output_typeP, buscoProteinOutput)
        if (copyToFtp) {
            COPY_PROTEIN_OUTPUT(buscoProteinSummaryOutput)
        }

    }
}





