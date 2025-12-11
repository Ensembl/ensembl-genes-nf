#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CREATE_SAMPLESHEET } from '../modules/create_samplesheet.nf'
include { REPORT_CONSENSUS } from '../modules/report_consensus.nf'


workflow TRANSLON_CONSENSUS {
    bed_files_ch = Channel
        .fromPath("${params.bed_results_dir}/*/*.bed")
        .map { bed_file ->
            def tool = bed_file.parent.name
            def sample_name = bed_file.baseName
            return tuple(sample_name, tool, bed_file)
        }
    
    // Group by sample name to collect all tools for each sample
    input_ch = bed_files_ch
        .groupTuple(by: 0)
        .map { sample_name, tools, bed_files ->
            // Create metadata
            def meta = [
                id: sample_name,
                tool_count: tools.size()
            ]
            
            // Sort tools and files together to maintain correspondence
            def sorted_pairs = [tools, bed_files].transpose().sort { it[0] }
            def sorted_tools = sorted_pairs.collect { it[0] }
            def sorted_files = sorted_pairs.collect { it[1] }
            
            return tuple(meta, sorted_tools, sorted_files)
        }

    CREATE_SAMPLESHEET(input_ch)

    REPORT_CONSENSUS(CREATE_SAMPLESHEET.out.results, params.ucsc_session_url)
}
