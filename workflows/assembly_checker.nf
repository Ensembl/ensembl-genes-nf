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

csvFile = file(params.csvFile)
  if( !csvFile.exists() ) {
    exit 1, "The specified csv file does not exist: ${params.csvfile}"
  }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.help) {
  log.info ''
  log.info 'Pipeline to check the quality of an assembly available on NCBI running Busco (auto-lineage)'
  log.info '-------------------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info '  nextflow -C ensembl-genes-nf/nextflow.config run ensembl-genes-nf/workflows/assembly_checker.nf --csvFile -profile slurm, standard'
  log.info ''
  log.info 'Options:'
  log.info '  --outDir STR              Output directory '
  log.info '  --csvFile STR             Path for the csv containing the db name'
  exit 1
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// MODULE: Loaded from modules/

include { PROCESS_ASSEMBLY } from '../modules/process_assembly.nf'
include { BUSCO_GENOME_AUTOLINEAGE } from '../modules/busco_genome_autolineage.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
        
        csvData = Channel.fromPath(params.csvFile).splitCsv(sep:',')
        
        //
        // MODULE: Download genomic sequences from NCBI FTP and store in the assembly accession directory
        //        
        
        PROCESS_ASSEMBLY (csvData)

        //
        // MODULE: Run Busco in genome mode and store the result in the assembly accession directory
        //
        BUSCO_GENOME_AUTOLINEAGE (PROCESS_ASSEMBLY.out.gca, PROCESS_ASSEMBLY.out.genome_file.flatten())
}
