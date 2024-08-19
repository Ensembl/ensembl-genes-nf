// main.nf
include {FETCH_GENOME} from "../modules/fetch_genome.nf"
include {CHECK_AND_DOWNLOAD_REPEATMODELER} from "../modules/check_and_download_rmlibrary.nf"
include {GENERATE_REPEATMODELER_LIBRARY} from "../modules/generate_repeatmodeler_library.nf"
inlcude {RUN_REPEATMASKER} from "../modules/run_repeatmasker.nf"

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

    // Run the FETCH_GENOME process
    data_with_gca
    .map { row -> row.gca }
    .set { gca_ch }

    FETCH_GENOME(gca_ch)

    // Create URLs for repeatmodeler files,check if they exist, and download them
    data_with_gca
    .map { row ->
        def url = "${params.repeats_ftp_base}/${row.species_name}/${row.gca}.repeatmodeler.fa"
        return [url, row.gca]
    }
    .set { url_ch }

   // Check and download RepeatModeler library
    check_and_download_result = CHECK_AND_DOWNLOAD_RMLIBRARY(url_ch)

    // Generate RepeatModeler library if it does not exist
    generate_library_result = check_and_download_result
        .filter { it.exists == false }
        .map { it.gca }
        .set { missing_library_gca_ch }

    GENERATE_REPEATMODELER_LIBRARY(missing_library_gca_ch)

    // Run RepeatMasker if the library exists or after generating the library
    repeatmasker_input = check_and_download_result
        .filter { it.exists == true }
        .map { it.gca }
        .merge(generate_library_result.map { it.gca })
        .set { repeatmasker_gca_ch }

    RUN_REPEATMASKER(repeatmasker_gca_ch)
}
