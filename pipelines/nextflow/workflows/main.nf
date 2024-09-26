nextflow.enable.dsl=2

include {FETCH_GENOME} from '../modules/fetch_genome.nf'
include {CHECK_AND_DOWNLOAD_RMLIBRARY} from '../modules/check_and_download_rmlibrary.nf'
include {GENERATE_REPEATMODELER_LIBRARY} from '../modules/generate_repeatmodeler_library.nf'
include {RUN_REPEATMASKER} from '../modules/run_repeatmasker.nf'

workflow {
    // Check if outDir parameter is defined
    if (!params.outDir) {
        error "Undefined --outDir parameter. Please provide the output directory's path"
    }

    // Check if csvFile parameter is defined and exists
    if (params.csvFile) {
        try {
            csvFile = file(params.csvFile, checkIfExists: true)
            println("CSV File exists: ${csvFile}")
        } catch (Exception e) {
            error "CSV file not found or could not be read: ${params.csvFile}"
        }
    } else {
        error 'CSV file not specified!'
    }

    // Attempt to read the CSV file
    data = Channel.fromPath(csvFile, type: 'file')
           .splitCsv(sep: ',', header: false)
           .map { row ->
                // Check for row correctness
                if (!row || row.size() < 2) {
                    error "Empty or malformed row in CSV: $row"
                }
                println("Processing CSV row: $row")  // Debug print each row of CSV
                return [species_name: row[0], gca: row[1]]
           }
           .set { data_with_gca }

    // Debugging: Print all rows of data_with_gca
    data_with_gca.view { row ->
        println("Data with GCA: $row")
    }

    // Check if channel data_with_gca is properly set
    if (!data_with_gca) {
        error "Failed to initialize data_with_gca channel."
    } else {
        println("data_with_gca channel initialized successfully.")
    }

    // Continue with the workflow logic below as previously
    data_with_gca
    .map { row -> row.gca }
    .set { gca_ch }

    FETCH_GENOME(gca_ch)

    data_with_gca
    .map { row ->
        def url = "${params.repeats_ftp_base}/${row.species_name}/${row.gca}.repeatmodeler.fa"
        println("Generated URL: $url")
        return [url, row.gca]
    }
    .set { url_ch }

    // Check if url_ch is set correctly
    url_ch.view { url ->
        println("URL Channel Content: $url")
    }

    // Check and download RepeatModeler library
    check_and_download_result = CHECK_AND_DOWNLOAD_RMLIBRARY(url_ch)
    check_and_download_result.view { result ->
        println("Check and Download Result: $result")
    }

    // Continue the workflow
    generate_library_result = check_and_download_result
        .filter { it.exists == false }
        .map { it.gca }
        .set { missing_library_gca_ch }

    GENERATE_REPEATMODELER_LIBRARY(missing_library_gca_ch)

    repeatmasker_input = check_and_download_result
        .filter { it.exists == true }
        .map { it.gca }
        .merge(generate_library_result.map { it.gca })
        .set { repeatmasker_gca_ch }

    RUN_REPEATMASKER(repeatmasker_gca_ch)
}
