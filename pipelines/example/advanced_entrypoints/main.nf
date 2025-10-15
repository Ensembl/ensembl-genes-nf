#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RUN_TOOLS } from './subworkflows_run_tools'
include { COMBINE_AND_COUNT } from './subworkflows_combine_and_count'

// Import the dynamic entry point system
import EntryPoints

// Parameters
params.outdir = 'test_results'
params.entry_point = 'auto'  // 'auto', 'FULL', 'RUN_TOOLS', 'COMBINE_AND_COUNT'
params.help = false

// Show help
if (params.help) {
    log.info """
    ${EntryPoints.getEntryPointsInfo()}

    Parameters:
      --entry_point    Entry point to start from (default: 'auto')
                       Options: auto, FULL, RUN_TOOLS, COMBINE_AND_COUNT
      --outdir         Output directory (default: 'test_results')
      --help           Show this help message

    Examples:
      # Automatic detection (recommended)
      nextflow run main.nf --outdir results

      # Force full workflow
      nextflow run main.nf --entry_point FULL --outdir results

      # Start from combine (requires tool outputs in --outdir)
      nextflow run main.nf --entry_point COMBINE_AND_COUNT --outdir results
    """
    System.exit(0)
}

// Helper to create channel from published files
def createChannelFromPublished(String baseDir, List<Map> requirements) {
    def allFiles = []

    requirements.each { req ->
        def pattern = "${baseDir}/${req.name}/${req.pattern}"
        def toolFiles = file(pattern)
        if (toolFiles.size() > 0) {
            log.info "  * Found ${toolFiles.size()} files in ${req.name}/"
            allFiles.addAll(toolFiles)
        }
    }

    return Channel.fromPath(allFiles)
        .map { file ->
            // Extract sample ID from filename
            def fileName = file.name
            def sampleId = fileName.replaceAll(/_(A|B|C|D|combined|line_count)\.txt$/, '')
            [[id: sampleId], file]
        }
}

// Main workflow with dynamic entry point resolution
workflow {
    // Determine which entry point to use
    def entryPoint = params.entry_point
    if (entryPoint == 'auto') {
        entryPoint = EntryPoints.determineEntryPoint(params.outdir)
        log.info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        log.info "Auto-detected entry point: ${entryPoint}"
        log.info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    } else {
        log.info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        log.info "User-specified entry point: ${entryPoint}"
        log.info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    }

    // Validate the entry point
    def validation = EntryPoints.validateEntryPoint(params.outdir, entryPoint)

    if (!validation.valid) {
        log.error "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        log.error "Entry Point Validation Failed"
        log.error "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        log.error validation.error
        log.error "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        error(validation.error)
    }

    log.info validation.message
    log.info ""

    // Create input channel
    input_ch = Channel.of(
        [id: 'sample1'],
        [id: 'sample2']
    )
    .map { meta ->
        // Create a dummy input file
        def dummy_file = file("${workflow.workDir}/dummy_input.txt")
        dummy_file.text = "dummy input"
        [meta, dummy_file]
    }

    // Execute the appropriate workflow based on entry point
    switch (entryPoint) {
        case 'FULL':
            log.info "Running complete workflow..."
            RUN_TOOLS(input_ch)
            COMBINE_AND_COUNT(RUN_TOOLS.out.results)
            break

        case 'RUN_TOOLS':
            log.info "Running tools only..."
            RUN_TOOLS(input_ch)
            break

        case 'COMBINE_AND_COUNT':
            log.info "Loading published files..."
            def requirements = EntryPoints.SUBWORKFLOWS[entryPoint].requires
            results_ch = createChannelFromPublished(params.outdir, requirements)
            log.info ""
            log.info "Running combine and count..."
            COMBINE_AND_COUNT(results_ch)
            break

        default:
            error "Unknown entry point: ${entryPoint}"
    }
}
