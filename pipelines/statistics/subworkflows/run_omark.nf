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

nextflow.enable.dsl = 2


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RUN OMARK WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In this subworkflow we fetch core database metadata, fetch protein sequences from the core database,
run Omamer to get orthologous groups, run OMark to get metakeys, and populate the core database with the results.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DB_METADATA } from '../modules/db_metadata.nf'
include { FETCH_PROTEINS } from '../modules/fetch_proteins.nf'
include { OMAMER_HOG } from '../modules/omamer_hog.nf'
include { OMARK } from '../modules/omark.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN OMARK WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RUN_OMARK {
    take:
    csvFile

    main:
    def ch_versions_file = channel.empty()
    // Read data from the CSV file, split it, and map each row to extract GCA and taxon values
    def data = channel.fromPath(csvFile, type: 'file', checkIfExists: true)
        .splitCsv(sep: ',', header: true)
        .map { row ->
            [
                gca: 'UNKNOWN',
                taxon_id: 'UNKNOWN',
                dbname: row.get('dbname'),
                species_id: row.get('species_id') ? row.get('species_id') : 1,
                protein_file: row.get('protein_file'),
            ]
        }
    def metadata = DB_METADATA(data).metadata.map { meta, metadata_file ->
        def lines = metadata_file.text.readLines()
        def new_taxon = lines[0].split('=')[1]
        def new_gca = lines[1].split('=')[1]
        def production_name = lines[2].split('=')[1]
        def updated_meta = meta + [taxon_id: new_taxon, gca: new_gca, production_species: production_name]
        updated_meta
    }
    ch_versions_file = ch_versions_file.mix(DB_METADATA.out.versions_file)
    // MODULE: Get canonical protein from db
    // 
    //def proteinData = FETCH_PROTEINS (db_meta).output.protein_file_output
    def proteinData = FETCH_PROTEINS(metadata).protein_file_output
    ch_versions_file = ch_versions_file.mix(FETCH_PROTEINS.out.versions_file)
    //
    // MODULE: Get orthologous groups from Omamer db 
    //
    def omamerData = OMAMER_HOG(proteinData).omamer_hog_output
    ch_versions_file = ch_versions_file.mix(OMAMER_HOG.out.versions_file)
    //
    // MODULE: Run Omark
    //        
    OMARK(omamerData)
    ch_versions_file = ch_versions_file.mix(OMARK.out.versions_file)

    emit:
    versions = ch_versions_file
}
