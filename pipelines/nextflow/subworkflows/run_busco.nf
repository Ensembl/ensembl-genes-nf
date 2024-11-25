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

include { BUSCO_CORE_METAKEYS as BUSCO_CORE_METAKEYS_PROTEIN } from '../modules/busco/busco_core_metakeys.nf'
include { BUSCO_CORE_METAKEYS as BUSCO_CORE_METAKEYS_GENOME } from '../modules/busco/busco_core_metakeys.nf'
include { BUSCO_DATASET } from '../modules/busco/busco_dataset.nf'
include { BUSCO_LINEAGES as BUSCO_GENOME_LINEAGE } from '../modules/busco/busco_lineages.nf'
include { BUSCO_LINEAGES as BUSCO_PROTEIN_LINEAGE } from '../modules/busco/busco_lineages.nf'
include { BUSCO_OUTPUT as BUSCO_GENOME_OUTPUT } from '../modules/busco/busco_output.nf'
include { BUSCO_OUTPUT as BUSCO_PROTEIN_OUTPUT } from '../modules/busco/busco_output.nf'
include { COPY_OUTPUT_TO_ENSEMBL_FTP as COPY_GENOME_OUTPUT } from '../modules/copy_output_to_ensembl_ftp.nf'
include { COPY_OUTPUT_TO_ENSEMBL_FTP as COPY_PROTEIN_OUTPUT } from '../modules/copy_output_to_ensembl_ftp.nf'
include { FETCH_GENOME } from '../modules/fetch_genome.nf'
include { FETCH_PROTEINS } from '../modules/fetch_proteins.nf'

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
        orthdb_clade_set = BUSCO_DATASET(db_meta).clade_dataset

        // // // Run Busco in genome mode
        if (busco_mode.contains('genome')) {
            def output_typeG = "genome"

            // Download genome via ncbi dataset API 
            FETCH_GENOME (db_meta)
            genome_fasta = FETCH_GENOME.out.genome_fasta

            buscoGenomeOutput = BUSCO_GENOME_LINEAGE(db_meta, orthdb_clade_set, output_typeG, genome_fasta)

            buscoGenomeSummaryOutput = BUSCO_GENOME_OUTPUT(db_meta, output_typeG, buscoGenomeOutput)

            // Copy BUSCO summary stats to ensembl FTP
            if (params.copyToFtp) {
                COPY_GENOME_OUTPUT(buscoGenomeSummaryOutput)
            }

            // Make and apply busco summary meta_keys patch directly to core:
            if(params.apply_busco_metakeys){
                BUSCO_CORE_METAKEYS_GENOME(buscoGenomeSummaryOutput)

            }
        }

        // // Run Busco in protein mode
        if (busco_mode.contains('protein')) {

            output_typeP = "protein"

            // Dump protein translations
            FETCH_PROTEINS (db_meta)
            protein_translations = FETCH_PROTEINS.out.translation_seqs

            buscoProteinOutput = BUSCO_PROTEIN_LINEAGE(db_meta, orthdb_clade_set, output_typeP, protein_translations)\

            buscoProteinSummaryOutput = BUSCO_PROTEIN_OUTPUT(db_meta, output_typeP, buscoProteinOutput)

            // Copy BUSCO summary stats to ensembl FTP
            if (params.copyToFtp) {
                COPY_PROTEIN_OUTPUT(buscoProteinSummaryOutput)
            }

            // Make and apply busco summary meta_keys patch directly to core:
            if(params.apply_busco_metakeys){
                BUSCO_CORE_METAKEYS_PROTEIN(buscoProteinSummaryOutput)
            }

        }
}



        // // // // Run Busco in genome mode
        // if (busco_mode.contains('genome')) {
        //     def output_typeG = "genome"

        //     // Download genome via ncbi dataset API 
        //     FETCH_GENOME(db_meta, orthdb_clade_set)
        //     expanded_genome_meta = FETCH_GENOME.out.meta
        //     genome_fasta = FETCH_GENOME.out.genome_fasta

        //     buscoGenomeOutput = BUSCO_GENOME_LINEAGE(expanded_genome_meta, output_typeG, genome_fasta)\
        //         .busco_report_output

        //     buscoGenomeSummaryOutput = BUSCO_GENOME_OUTPUT(output_typeG, expanded_genome_meta, buscoGenomeOutput)

        //     // Copy BUSCO summary stats to ensembl FTP
        //     if (params.copyToFtp) {
        //         COPY_GENOME_OUTPUT(buscoGenomeSummaryOutput)
        //     }

        //     // Make and apply busco summary meta_keys patch directly to core:
        //     if(params.apply_busco_metakeys){
        //         BUSCO_CORE_METAKEYS_GENOME(buscoGenomeSummaryOutput)

        //     }
        // }