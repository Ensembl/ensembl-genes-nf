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
if (!params.mode) {
  exit 1, "Undefined --mode parameter. Please define Busco running mode"
}

if (params.csvFile) {
    ch_csvFile = file(params.csvFile, checkIfExists: true)
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
  exit 1
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// MODULE: Loaded from modules/

include { BUSCO_DATASET } from '../modules/busco_dataset.nf'
include { SPECIES_OUTDIR } from '../modules/species_outdir.nf'
include { FETCH_GENOME } from '../modules/fetch_genome.nf'
include { FETCH_PROTEINS } from '../modules/fetch_proteins.nf'
include { BUSCO_GENOME_LINEAGE } from '../modules/busco_genome_lineage.nf'
include { BUSCO_PROTEIN_LINEAGE } from '../modules/busco_protein_lineage.nf'
include { BUSCO_GENOME_OUTPUT } from '../modules/busco_genome_output.nf'
include { BUSCO_PROTEIN_OUTPUT } from '../modules/busco_protein_output.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
        csvData = ch_csvFile.splitCsv()
        buscoModes = Channel.fromList(busco_mode)

        //
        // MODULE: Get the closest Busco dataset from the taxonomy classification stored in db meta table 
        //        
        BUSCO_DATASET (csvData.flatten())
        
        //
        // MODULE: Create directory path for FTP
        //
        SPECIES_OUTDIR (BUSCO_DATASET.out.dbname, BUSCO_DATASET.out.busco_dataset)
        SPECIES_OUTDIR.out.combine(buscoModes).branch {
                        protein: it[3] == 'protein'
                        genome: it[3] == 'genome'
             }.set { ch_mode }
        
        //
        // MODULE: Get genomic sequences from db
        //        
        FETCH_GENOME (ch_mode.genome)
        
        //
        // MODULE: Run Busco in genome mode
        //        
        BUSCO_GENOME_LINEAGE (FETCH_GENOME.out.fasta.flatten(), FETCH_GENOME.out.output_dir, FETCH_GENOME.out.db_name, FETCH_GENOME.out.busco_dataset)

        //
        // MODULE: Edit Busco summary file
        //
        BUSCO_GENOME_OUTPUT(BUSCO_GENOME_LINEAGE.out.species_outdir)        
        
        //
        // MODULE: Get canonical protein from db
        //        
        FETCH_PROTEINS (ch_mode.protein)
        
        //
        // MODULE: Run Busco in protein mode
        //        
        BUSCO_PROTEIN_LINEAGE (FETCH_PROTEINS.out.fasta.flatten(), FETCH_PROTEINS.out.output_dir, FETCH_PROTEINS.out.db_name, FETCH_PROTEINS.out.busco_dataset)
        
        //
        // MODULE: Edit Busco summary file
        //        
        BUSCO_PROTEIN_OUTPUT(BUSCO_PROTEIN_LINEAGE.out.species_outdir) 
}
