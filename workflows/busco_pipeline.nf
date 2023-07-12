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
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   VALIDATE INPUTS
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (!params.host) {
  exit 1, "Undefined --host parameter. Please provide the server host for the db connection"
}

if (!params.port) {
  exit 1, "Undefined --port parameter. Please provide the server port for the db connection"
}
if (!params.user) {
  exit 1, "Undefined --user parameter. Please provide the server user for the db connection"
}

if (!params.enscode) {
  exit 1, "Undefined --enscode parameter. Please provide the enscode path"
}
if (!params.outDir) {
  exit 1, "Undefined --outDir parameter. Please provide the output directory's path"
}
if (!params.cacheDir) {
  exit 1, "Undefined --cacheDir parameter. Please provide the cache dir directory's path"
}
if (!params.mode) {
  exit 1, "Undefined --mode parameter. Please define Busco running mode"
}

if (params.csvFile) {
    csvFile = file(params.csvFile, checkIfExists: true)
} else {
    exit 1, 'CSV file not specified!'
}

busco_mode = []
if (params.mode instanceof java.lang.String) {
  busco_mode = [params.mode]
}
else {
  busco_mode = params.mode
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.help) {
  log.info ''
  log.info 'Pipeline to run Busco score in protein and/or genome mode'
  log.info '-------------------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info '  nextflow -C ensembl-genes-nf/nextflow.config run ensembl-genes-nf/iworkflows/busco_pipeline.nf --enscode --csvFile --mode'
  log.info ''
  log.info 'Options:'
  log.info '  --host STR                   Db host server '
  log.info '  --port INT                   Db port  '
  log.info '  --user STR                   Db user  '
  log.info '  --enscode STR                Enscode path '
  log.info '  --outDir STR                 Output directory. Default is workDir'
  log.info '  --csvFile STR                Path for the csv containing the db name'
  log.info '  --mode STR                   Busco mode: genome or protein, default is to run both'
  log.info '  --bioperl STR                BioPerl path (optional)'
  log.info '  --project STR                Project, for the formatting of the output ("ensembl" or "brc")'
  exit 1
}

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
include { BUSCO_OUTPUT as BUSCO_GENOME_OUTPUT } from '../modules/busco_output.nf'
include { BUSCO_OUTPUT as BUSCO_PROTEIN_OUTPUT } from '../modules/busco_output.nf'
include { FASTA_OUTPUT as FASTA_GENOME_OUTPUT } from '../modules/fasta_output.nf'
include { FASTA_OUTPUT as FASTA_PROTEIN_OUTPUT } from '../modules/fasta_output.nf'
include { SPECIES_METADATA } from '../modules/species_metadata.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    csvData = Channel.fromPath(params.csvFile).splitCsv()

    // Get db name and its metadata
    db = csvData.flatten()
    db_meta = SPECIES_METADATA(db, params.outDir, params.project)
      .splitCsv(header: true)

    // Get the closest Busco dataset from the taxonomy classification stored in db meta table 
    db_dataset = BUSCO_DATASET(db_meta)
    
    // Run Busco in genome mode
    if (busco_mode.contains('genome')) {
        genome_data = FETCH_GENOME(db_dataset, params.cacheDir)
        busco_genome_output = BUSCO_GENOME_LINEAGE(genome_data)
        BUSCO_GENOME_OUTPUT(busco_genome_output, "genome", params.project)
        if (params.project == 'ensembl') {
          FASTA_GENOME_OUTPUT(genome_data, params.project, 'genome')
        }
    }
    
    // Run Busco in protein mode
    if (busco_mode.contains('protein')) {
        protein_data = FETCH_PROTEINS (db_dataset, params.cacheDir)
        busco_protein_output = BUSCO_PROTEIN_LINEAGE(protein_data)
        BUSCO_PROTEIN_OUTPUT(busco_protein_output, "protein", params.project)
        if (params.project == 'ensembl') {
          FASTA_PROTEIN_OUTPUT(protein_data, params.project, 'fasta')
        }
    }
}
