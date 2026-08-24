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
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DB_METADATA } from '../modules/db_metadata.nf'
include { BUSCO_DATASET } from '../modules/busco_dataset.nf'
include { FETCH_GENOME } from '../modules/fetch_genome.nf'
include { FETCH_PROTEINS } from '../modules/fetch_proteins.nf'
include { BUSCO_LINEAGE } from '../modules/busco_lineage.nf'
include { BUSCO_CORE_METAKEYS } from '../modules/busco_core_metakeys.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN BUSCO SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In this subworkflow we run BUSCO on either genome or protein mode based on user input.
We first fetch metadata from the core database or from NCBI depending on the input parameters.
We then select the appropriate BUSCO dataset based on the taxon_id.
Finally, we run BUSCO in the selected mode(s) and populate the core database with the results.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
*/
workflow RUN_BUSCO {
    take:
    csvFile
    protein_input
    use_shared_proteins

    main:
    def data
    def busco_mode = params.busco_mode == 'both' ? ['protein', 'genome'] : [params.busco_mode]
    if (params.run_busco_ncbi && !params.run_busco_core) {
        busco_mode = ['genome']
        data = channel.fromPath(csvFile, type: 'file', checkIfExists: true)
            .splitCsv(sep: ',', header: true)
            .map { row ->
                [
                    gca: row.get('gca'),
                    taxon_id: row.get('taxon_id'),
                    dbname: 'UNKNOWN',
                    species_id: row.get('species_id') ? row.get('species_id') : 1,
                    busco_mode: busco_mode,
                    busco_dataset: row.get('busco_dataset'),
                    genome_file: row.get('genome_file'),
                    protein_file: row.get('protein_file'),
                ]
            }
    }
    else if (params.run_busco_core) {
        data = channel.fromPath(csvFile, type: 'file', checkIfExists: true)
            .splitCsv(sep: ',', header: true)
            .map { row ->
                [
                    gca: 'UNKNOWN',
                    taxon_id: 'UNKNOWN',
                    dbname: row.get('dbname'),
                    species_id: row.get('species_id') ? row.get('species_id') : 1,
                    busco_mode: busco_mode,
                    busco_dataset: row.get('busco_dataset'),
                    genome_file: row.get('genome_file'),
                    protein_file: row.get('protein_file'),
                ]
            }
    }
    else {
        error("At least one of the parameters params.run_busco_ncbi or params.run_busco_core must be set to true")
    }

    data.view { d -> "GCA: ${d.gca}, Taxon ID: ${d.taxon_id}, Core name: ${d.dbname}, Species ID: ${d.species_id}" }
    // Get metadata from the database
    def metadata = DB_METADATA(data).metadata.map { meta, metadata_file ->
        def lines = metadata_file.text.readLines()
        def new_taxon = lines[0].split('=')[1]
        def new_gca = lines[1].split('=')[1]
        def production_name = lines[2].split('=')[1]
        def updated_meta = meta + [taxon_id: new_taxon, gca: new_gca, production_name: production_name]
        updated_meta
    }
    ch_versions_file = channel.empty()
    ch_versions_file = ch_versions_file.mix(DB_METADATA.out.versions_file)
    // Get the closest Busco dataset from the taxonomy classification stored in db meta table

    def dataset_db = BUSCO_DATASET(metadata).busco_dataset_output
        .view { item -> "Channel contains: ${item}" }
        .map { tuple_meta, stdout_file ->
            def new_meta = tuple_meta + [busco_dataset: tuple_meta.busco_dataset ? tuple_meta.busco_dataset.trim() : stdout_file.trim()]
            return new_meta
        }
        .view { meta -> "Mapped result: gca=${meta.gca}, dbname=${meta.dbname}, busco_dataset=${meta.busco_dataset}" }
    ch_versions_file = ch_versions_file.mix(BUSCO_DATASET.out.versions_file)

    def buscoInput
    if (busco_mode.contains('genome')) {
        def genomeData = FETCH_GENOME(dataset_db).fasta_file_output
            .map { meta, fasta_file ->
                return tuple(meta + [busco_mode: 'genome'], fasta_file)
            }
            .view { meta, fasta_file -> "Mapped result: gca=${meta.gca}, dbname=${meta.dbname}, busco_dataset=${meta.busco_dataset}, file=${fasta_file}" }
        ch_versions_file = ch_versions_file.mix(FETCH_GENOME.out.versions_file)
        buscoInput = genomeData
    }

    if (busco_mode.contains('protein')) {
        def proteinData
        if (use_shared_proteins) {
            proteinData = protein_input
                .map { meta, fasta_file -> tuple(meta.dbname, fasta_file) }
                .join(dataset_db.map { meta -> tuple(meta.dbname, meta) })
                .map { row -> tuple(row[2], row[1]) }
        } else {
            proteinData = FETCH_PROTEINS(dataset_db).fasta_file_output
        }
        proteinData = proteinData.map { meta, fasta_file ->
                return tuple(meta + [busco_mode: 'protein'], fasta_file)
        }
        if (!use_shared_proteins) {
            ch_versions_file = ch_versions_file.mix(FETCH_PROTEINS.out.versions_file)
        }

        buscoInput = buscoInput != null ? buscoInput.mix(proteinData) : proteinData
    }

    if (buscoInput == null) {
        error("Invalid BUSCO mode: ${params.busco_mode}")
    }

    def buscoOutput = BUSCO_LINEAGE(buscoInput).busco_lineage_output
    ch_versions_file = ch_versions_file.mix(BUSCO_LINEAGE.out.versions_file)

    BUSCO_CORE_METAKEYS(buscoOutput)
    ch_versions_file = ch_versions_file.mix(BUSCO_CORE_METAKEYS.out.versions_file)

    emit:
    versions = ch_versions_file
}
