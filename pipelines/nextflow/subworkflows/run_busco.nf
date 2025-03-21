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
// include { COPY_OUTPUT_TO_ENSEMBL_FTP as COPY_GENOME_OUTPUT } from '../modules/copy_output_to_ensembl_ftp.nf'
// include { COPY_OUTPUT_TO_ENSEMBL_FTP as COPY_PROTEIN_OUTPUT } from '../modules/copy_output_to_ensembl_ftp.nf'
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
        // Get the closest BUSCO dataset from the taxonomy classification stored in db meta table
        BUSCO_DATASET(db_meta).clade_dataset
            .map { row -> [
                insdc_acc: row[0], taxonomy_id: row[1], core: row[2], production_name: row[3],
                organism_name: row[4], annotation_source: row[5], ortho_db: row[6]
            ]}
            .set{ orthodb_amended_meta }

        // Run BUSCO in genome mode
        if (busco_mode.contains('genome')) {

            output_typeG = "genome"

            // Download genome via ncbi dataset API
            FETCH_GENOME (orthodb_amended_meta).genome_fasta
                .map { row -> [
                    insdc_acc: row[0], taxonomy_id: row[1], core: row[2], production_name: row[3],
                    organism_name: row[4], annotation_source: row[5], ortho_db: row[6], genome: row[7]
                ]}
                .set{ metadata_genome_fna }

            BUSCO_GENOME_LINEAGE(metadata_genome_fna, output_typeG).busco_report_output
                .map { row -> [
                    insdc_acc: row[0], taxonomy_id: row[1], core: row[2],production_name: row[3],
                    organism_name: row[4], annotation_source: row[5], report: row[6]
                ]}
                .set{ buscoGenomeOutput }

            buscoGenomeSummaryOutput = BUSCO_GENOME_OUTPUT(buscoGenomeOutput, output_typeG)

            // Copy BUSCO summary stats to ensembl FTP
            // if (copyToFtp) {
            //     COPY_GENOME_OUTPUT(buscoGenomeSummaryOutput)
            // }

            // Make and apply busco summary meta_keys patch directly to core:
            if(params.apply_busco_metakeys){
                BUSCO_CORE_METAKEYS_GENOME(buscoGenomeSummaryOutput)
            }
        }

        // Run Busco in protein mode
        if (busco_mode.contains('protein')) {

            output_typeP = "protein"

            // Dump protein translations
            FETCH_PROTEINS (orthodb_amended_meta).translations
                .map { row -> [
                    insdc_acc: row[0], taxonomy_id: row[1], core: row[2], production_name: row[3],
                    organism_name: row[4], annotation_source: row[5], ortho_db: row[6], translations: row[7]
                ]}
                .set{ metadata_prot_trans }

            BUSCO_PROTEIN_LINEAGE(metadata_prot_trans, output_typeP).busco_report_output
                .map { row -> [
                    insdc_acc: row[0], taxonomy_id: row[1], core: row[2], production_name: row[3],
                    organism_name: row[4], annotation_source: row[5], report: row[6]
                ]}
                .set{ buscoProteinOutput }

            buscoProteinSummaryOutput = BUSCO_PROTEIN_OUTPUT(buscoProteinOutput, output_typeP)

            // Copy BUSCO summary stats to ensembl FTP
            // if (copyToFtp) {
            //     COPY_PROTEIN_OUTPUT(buscoProteinSummaryOutput)
            // }

            // Make and apply busco summary meta_keys patch directly to core:
            if(params.apply_busco_metakeys){
                BUSCO_CORE_METAKEYS_PROTEIN(buscoProteinSummaryOutput)
            }

        }
}