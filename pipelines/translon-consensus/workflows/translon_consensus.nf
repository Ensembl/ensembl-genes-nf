#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RENAME_BED } from '../modules/rename_bed.nf'
include { CREATE_SAMPLESHEET } from '../modules/create_samplesheet.nf'
include { REPORT_CONSENSUS } from '../modules/report_consensus.nf'


workflow TRANSLON_CONSENSUS {
    bed_files_ch = Channel
        .fromPath("${params.bed_results_dir}/*/*.bed")
        .map { bed_file ->
            def tool = bed_file.parent.name
            def sample_name = bed_file.baseName
            def meta = [id: sample_name]
            return tuple(meta, tool, bed_file)
        }
    
    // Rename files to include tool name
    RENAME_BED(bed_files_ch)
    
    // Group by sample name to collect all renamed files
    input_ch = RENAME_BED.out
        .groupTuple(by: 0)
        .map { meta, renamed_files ->
            // Sort files by name for consistency
            def sorted_files = renamed_files.sort { it.name }
            return tuple(meta, sorted_files)
        }

    CREATE_SAMPLESHEET(input_ch)

    REPORT_CONSENSUS(CREATE_SAMPLESHEET.out, params.ucsc_session_url)
}
