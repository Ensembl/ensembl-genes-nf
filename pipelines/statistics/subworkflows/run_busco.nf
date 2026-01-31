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

include { DB_METADATA } from '../modules/db_metadata.nf'
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
    def data 
    def busco_mode = params.busco_mode == 'both' ? ['protein', 'genome'] : [params.busco_mode]
    if(params.run_busco_ncbi && !params.run_busco_core){
        busco_mode = 'genome'
        data = channel.fromPath(csvFile, type: 'file', checkIfExists: true)
                .splitCsv(sep:',', header:true)
                .map { row -> 
                    [gca:row.get('gca'), 
                    taxon_id:row.get('taxon_id'), 
                    dbname:'UNKNOWN', 
                    species_id:row.get('species_id') ? row.get('species_id'):1,
                    busco_mode:busco_mode,
                    busco_dataset:row.get('busco_dataset'),
                    genome_file:row.get('genome_file'),
                    protein_file:row.get('protein_file')
                    ]
                    }
    }
    else if(params.run_busco_core){
        data = channel.fromPath(params.csvFile, type: 'file', checkIfExists: true)
                .splitCsv(sep:',', header:true)
                .map { row -> [
                    gca:'UNKNOWN', 
                    taxon_id:'UNKNOWN', 
                    dbname:row.get('dbname'), 
                    species_id:row.get('species_id') ? row.get('species_id'):1,
                    busco_mode:busco_mode,
                    busco_dataset:row.get('busco_dataset'),
                    genome_file:row.get('genome_file'),
                    protein_file:row.get('protein_file')]}
    }
    else{
        error "At least one of the parameters params.run_busco_ncbi or params.run_busco_core must be set to true"
    }
    
    data.view { d -> "GCA: ${d.gca}, Taxon ID: ${d.taxon_id}, Core name: ${d.dbname}, Species ID: ${d.species_id}" }
    // Get metadata from the database
    def metadata = DB_METADATA(data).metadata
    .map { meta, metadata_file ->
            def lines = metadata_file.text.readLines()
            def new_taxon = lines[0].split('=')[1]
            def new_gca = lines[1].split('=')[1]
            def production_name = lines[2].split('=')[1]
            def updated_meta = meta + [taxon_id: new_taxon, gca: new_gca, production_species: production_name]
            updated_meta
        }
    ch_versions_file = channel.empty()
    ch_versions_file = ch_versions_file.mix(DB_METADATA.out.versions_file)
    // Get the closest Busco dataset from the taxonomy classification stored in db meta table

    def dataset_db = BUSCO_DATASET(metadata).busco_dataset_output
        .view { item -> "Channel contains: ${item}" }
        .map{ tuple_meta, stdout_file  ->
            def new_meta=tuple_meta + [ busco_dataset: tuple_meta.busco_dataset ? tuple_meta.busco_dataset.trim() : stdout_file.trim()]
            return new_meta
        }.view { meta -> "Mapped result: gca=${meta.gca}, dbname=${meta.dbname}, busco_dataset=${meta.busco_dataset}" }
    ch_versions_file = ch_versions_file.mix(BUSCO_DATASET.out.versions_file)
    // Run Busco in genome mode
    if (busco_mode.contains('genome')) {
        def genomeData = FETCH_GENOME(dataset_db).genome_file_output
        .map { meta, fna_file ->
            return tuple(meta, fna_file)
        }.view { meta,fna_file -> "Mapped result: gca=${meta.gca}, dbname=${meta.dbname}, busco_dataset=${meta.busco_dataset}, file=${fna_file}" }
        ch_versions_file = ch_versions_file.mix(FETCH_GENOME.out.versions_file)
        def buscoGenomeOutput = BUSCO_GENOME_LINEAGE(genomeData).busco_genome_lineage_output
        ch_versions_file = ch_versions_file.mix(BUSCO_GENOME_LINEAGE.out.versions_file)
        BUSCO_CORE_METAKEYS_GENOME(buscoGenomeOutput)
        ch_versions_file = ch_versions_file.mix(BUSCO_CORE_METAKEYS_GENOME.out.versions_file)
    }
    
    // Run Busco in protein mode
    if (busco_mode.contains('protein')) {
        def proteinData = FETCH_PROTEINS(dataset_db).protein_file_output
        ch_versions_file = ch_versions_file.mix(FETCH_PROTEINS.out.versions_file)
        def buscoProteinOutput = BUSCO_PROTEIN_LINEAGE(proteinData).busco_protein_lineage_output
        ch_versions_file = ch_versions_file.mix(BUSCO_PROTEIN_LINEAGE.out.versions_file)
        BUSCO_CORE_METAKEYS_PROTEIN(buscoProteinOutput)
        ch_versions_file = ch_versions_file.mix(BUSCO_CORE_METAKEYS_PROTEIN.out.versions_file)
    }
    emit:
    versions = ch_versions_file
}





