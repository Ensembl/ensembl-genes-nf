#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// CLI parameters
params.input   = null
params.db_pass = null

// Require input and db_pass
if (!params.input)   error "Please provide --input <csv_file>"
if (!params.db_pass) error "Please provide --db_pass <password>"

// Includes
include { BROKEN_TRANSLATIONS }      from './subworkflows/broken_translations.nf'
include { METADATA }                 from './subworkflows/metadata.nf'
include { MISSING_TRANSLATIONS }     from './subworkflows/missing_features'
include { FIX_GENE_DISPLAY_XREFS }   from './subworkflows/broken_xrefs.nf'



// Define workflow
workflow {
    // Read CSV -> Channel
    // Assuming CSV header: core_db,some_other_columns...
    core_info_ch = Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def core_db = row.database_name
            def gca = core_db.tokenize('_')[2]  // Extract the GCA from db name
            
            // Create meta map with id field (required by the process)
            def meta = [
                id: core_db,           // This is what the process uses in tag "${meta.id}"
                gca: gca,
                server: 'gb1-w'
            ]
            
            // Return tuple matching the expected input: tuple val(meta), val(core_name)
            tuple(meta, core_db)
        }
    // Run broken translations
    // broken_results = BROKEN_TRANSLATIONS(core_info_ch, params.db_pass)

    // Run metadata steps
    // METADATA(core_info_ch, 
    //          params.skip_dc) // Set to true to skip final verification datachecks

    // Run missing translations
    // missing_results = MISSING_TRANSLATIONS(core_info_ch, params.db_pass)

    FIX_GENE_DISPLAY_XREFS(core_info_ch)
}
