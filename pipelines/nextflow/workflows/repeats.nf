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
VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (!params.outDir) {
    exit 1, "Undefined --outDir parameter. Please provide the output directory's path"
}

if (params.csvFile) {
    csvFile = file(params.csvFile, checkIfExists: true)
} else {
    exit 1, 'CSV file not specified!'
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.help) {
    log.info ''
    log.info 'Pipeline to run RepeatModeler library generation'
    log.info '-------------------------------------------------------'
    log.info ''
    log.info 'Usage: '

    exit 1
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

inlcude { RUN_REPEATMASKER  } from '../modules/repeatmasker/repeatmasker.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow REPEATS{
        // Read data from the CSV file, split it, and map each row to extract GCA and species name values
        // fetch the genome from ncbi for the GCA
	data = Channel.fromPath(params.csvFile, type: 'file', checkIfExists: true)
	       .splitCsv(sep:',', header:true)
	       .map { row -> [gca:row.get('gca'), species_name:row.get('species_name')]
           .map { row ->
               def url = "${params.repeats_ftp_base}/${row.species_name}/${row.gca}.repeatmodeler.fa"
               return [row, url]
               }
               .set { data_with_url }
	def exists = CHECK_FILE_EXISTS(url)

//if the file repeatmodeler.fa exists downdload, run repeatmasker
    data_with_url.flatMap { row, url ->
        CHECK_FILE_EXISTS(url).map { exists ->
            def exists_file = exists.toFile().text.trim()
            if (exists_file == 'true') {
                return [url, row]
            }
            return null
        }
    }.filter { it != null }
    .set { existing_files }

    existing_files | FETCH_GENOME| DOWNLOAD_FILE | RUN_REPEATMASKER

    RUN_REPEATMASKER.out.view { result ->
        println("RepeatMasker output for ${result.row.gca}: ${result}")
    }
}
}
