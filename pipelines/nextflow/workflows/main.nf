nextflow.enable.dsl=2
import java.nio.file.Files

include {FETCH_GENOME} from '../modules/fetch_genome.nf'
include {FETCH_REPEAT_MODEL} from '../modules/check_and_download_rmlibrary.nf'
include {GENERATE_REPEATMODELER_LIBRARY} from '../modules/generate_repeatmodeler_library.nf'

workflow {
    // Check if outDir parameter is defined
    if (!params.outDir) {
        error "Undefined --outDir parameter. Please provide the output directory's path"
    }

    // Check if csvFile parameter is defined and exists
    if (params.csvFile) {
        csvFile = file(params.csvFile, checkIfExists: true)
    } else {
        error 'CSV file not specified!'
    }
// Read data from the CSV file, split it, and map each row to extract species_name and GCA values
    data = Channel.fromPath(csvFile, type: 'file')
           .splitCsv(sep: ',', header: false)
           .map { row -> [species_name: row[0], gca: row[1]] }
           .set { data_with_gca }
    FETCH_GENOME(data_with_gca)
        .set{fetched_genome}

   def fetch_repeat = FETCH_REPEAT_MODEL(fetched_genome)

// Filter tuples where repeatmodeler file is incomplete
incomplete_repeat_models = fetch_repeat.filter { species, gca, genome_file, repeatmodeler_file ->
    def fileContent = repeatmodeler_file.text
    return fileContent.contains("No repeatmodeler file available")
}

//Send the GCA and genome path to the GENERATE_REPEATMODELER_LIBRARY process
GENERATE_REPEATMODELER_LIBRARY(incomplete_repeat_models)
}    
