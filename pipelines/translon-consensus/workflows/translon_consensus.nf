#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RENAME_BED } from '../modules/rename_bed.nf'
include { CREATE_SAMPLESHEET } from '../modules/create_samplesheet.nf'
include { REPORT_CONSENSUS } from '../modules/report_consensus.nf'
include { GENERATE_HTML_REPORT } from '../modules/generate_report.nf'


workflow TRANSLON_CONSENSUS {
    bed_files_ch = Channel
        .fromPath("${params.bed_results_dir}/*/*.bed12")
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
            // Add tool count to metadata
            def enriched_meta = meta + [tool_count: sorted_files.size()]
            return tuple(enriched_meta, sorted_files)
        }

    // Split into single-tool and multi-tool samples
    multi_tool_ch = input_ch
        .filter { meta, files -> meta.tool_count >= 2 }

    single_tool_ch = input_ch
        .filter { meta, files -> meta.tool_count == 1 }

    CREATE_SAMPLESHEET(multi_tool_ch)

    // Prepare GTF input as value channel (use file or empty placeholder)
    // Value channel allows reuse across multiple samples
    gtf_ch = params.gencode_gtf ? Channel.value(file(params.gencode_gtf, checkIfExists: true)) : Channel.value(file('NO_FILE'))

    // Only run consensus analysis for samples with 2+ tools
    REPORT_CONSENSUS(CREATE_SAMPLESHEET.out, params.ucsc_session_url, gtf_ch)

    // Log single-tool samples (no consensus possible)
    single_tool_ch.subscribe { meta, files ->
        log.warn "Sample ${meta.id} has only 1 tool - skipping consensus analysis (requires 2+ tools)"
    }

    // Collect all consensus results and generate HTML report
    all_results = REPORT_CONSENSUS.out.results
        .map { meta, files -> files }
        .flatten()
        .collect()

    GENERATE_HTML_REPORT(all_results, params.ucsc_session_url)
}
