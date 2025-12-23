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

include { BUSCO_DATASET } from '../modules/busco_dataset.nf'
include { FETCH_GENOME } from '../modules/fetch_genome.nf'
include { FETCH_PROTEINS } from '../modules/fetch_proteins.nf'
include { BUSCO_GENOME_LINEAGE } from '../modules/busco_genome_lineage.nf'
include { BUSCO_PROTEIN_LINEAGE } from '../modules/busco_protein_lineage.nf'
//include { BUSCO_OUTPUT as BUSCO_GENOME_OUTPUT } from '../modules/busco_output.nf'
//include { BUSCO_OUTPUT as BUSCO_PROTEIN_OUTPUT } from '../modules/busco_output.nf'
include { BUSCO_CORE_METAKEYS as BUSCO_CORE_METAKEYS_PROTEIN } from '../modules/busco_core_metakeys.nf'
include { BUSCO_CORE_METAKEYS as BUSCO_CORE_METAKEYS_GENOME } from '../modules/busco_core_metakeys.nf'
//include { COPY_OUTPUT_TO_ENSEMBL_FTP as COPY_GENOME_OUTPUT } from '../modules/copy_output_to_ensembl_ftp.nf'
//include { COPY_OUTPUT_TO_ENSEMBL_FTP as COPY_PROTEIN_OUTPUT } from '../modules/copy_output_to_ensembl_ftp.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow RUN_BUSCO{
    take:                 
    csvFile

    main:
    //busco_script = file("${projectDir}/bin/busco_metakeys_patch.py")
    // Validate required parameters
    if (params.fetch == null) {
        error "params.fetch must be defined as true or false"
    }
    if(params.run_busco_ncbi && !params.run_busco_core){
        def busco_mode = 'genome'
        // Read data from the CSV file, split it, and map each row to extract GCA and taxon values
        def data = Channel.fromPath(csvFile, type: 'file', checkIfExists: true)
                .splitCsv(sep:',', header:true)
                .map { row -> 
                    [gca:row.get('gca'), 
                    taxon_id:row.get('taxon_id'), 
                    core:'UNKNOWN', 
                    busco_mode:busco_mode
                    ]
                    }
                
        
    }
    else if(params.run_busco_core){
        def busco_mode = params.busco_mode
        // Read data from the CSV file, split it, and map each row to extract GCA and taxon values
        data = Channel.fromPath(params.csvFile, type: 'file', checkIfExists: true)
                .splitCsv(sep:',', header:true)
                .map { row -> [
                    gca:'UNKNOWN', 
                    taxon_id:'UNKNOWN', 
                    core:row.get('core'), 
                    species_id:row.get('species_id'), busco_mode:busco_mode]}
                
        
    }
    else{
        error "At least one of the parameters params.run_busco_ncbi or params.run_busco_core must be set to true"
    }

    // Get the closest Busco dataset from the taxonomy classification stored in db meta table
    //def db_meta1=db_meta
    //db_meta1.flatten().view { d -> "GCA: ${d.gca}, Taxon ID: ${d.taxon_id}, Core name: ${d.core}, Species ID: ${d.species_id}" }
    //def buscoDataset = params.busco_dataset ? params.busco_dataset.trim() : meta.busco_dataset.trim() 

    def dataset_db = BUSCO_DATASET(data).output.busco_dataset_output.map{ tuple_meta, stdout ->
        [
        gca: tuple_meta.gca,
        core: tuple_meta.core,
        species_id: tuple_meta.species_id,
        busco_mode: tuple_meta.busco_mode,
        busco_dataset: params.busco_dataset ? params.busco_dataset.trim() :stdout.trim()
    ]
}

    // Run Busco in genome mode
    if (params.busco_mode.toLowerCase().contains('genome')) {
        //def output_typeG = "genome"
        def genomeData = FETCH_GENOME(dataset_db).output.genome_file_output
        def buscoGenomeOutput = BUSCO_GENOME_LINEAGE(genomeData).output.busco_genome_lineage_output
        BUSCO_CORE_METAKEYS_GENOME(buscoGenomeOutput)
        //if(params.apply_busco_metakeys){
            
        //}
    }
    
    // Run Busco in protein mode
    if (params.busco_mode.toLowerCase().contains('protein')) {
        //def output_typeP = "protein"
        def proteinData = FETCH_PROTEINS(dataset_db).output.protein_file_output
        def buscoProteinOutput = BUSCO_PROTEIN_LINEAGE(proteinData).output.busco_protein_lineage_output
        BUSCO_CORE_METAKEYS_PROTEIN(buscoProteinOutput)
        //def (buscoProteinSummaryOutput) = BUSCO_PROTEIN_OUTPUT(output_typeP, buscoProteinOutput)
        //if (copyToFtp) {
        //    COPY_PROTEIN_OUTPUT(buscoProteinSummaryOutput)
        //}
        //def buscoProteinSummaryOutput1=buscoProteinSummaryOutput
        //if(params.apply_busco_metakeys){
            
        //}

    }
}





