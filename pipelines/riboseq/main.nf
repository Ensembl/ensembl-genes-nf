#!/usr/bin/env nextflow
/*
========================================================================================
    RIBOSEQ PIPELINE
========================================================================================
    Ribosome profiling data processing pipeline
    - Data acquisition and read collapsing
    - Quality control
    - STAR alignment (genome + transcriptome)
    - RiboMetric and RiboWaltz analysis (offsets, QC, profiles)
    - BEDgraph and BigWig generation
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { validateParameters } from 'plugin/nf-schema'

// Validate parameters against schema
validateParameters()


/*
========================================================================================
    IMPORT SUBWORKFLOWS
========================================================================================
*/

include { DATA_ACQUISITION } from './subworkflows/data_acquisition.nf'
include { QUALITY_CONTROL } from './subworkflows/quality_control.nf'
include { ALIGNMENT } from './subworkflows/alignment.nf'
include { ANALYSIS } from './subworkflows/analysis.nf'
include { POST_PROCESSING } from './subworkflows/post_processing.nf'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    //
    // SUBWORKFLOW: Data acquisition and read collapsing
    //
    DATA_ACQUISITION(
        params.sample_sheet,
        file(params.adapter_list)
    )

    //
    // SUBWORKFLOW: Quality control of collapsed reads
    //
    QUALITY_CONTROL(
        DATA_ACQUISITION.out.samples
    )

    //
    // SUBWORKFLOW: STAR alignment (genome + transcriptome)
    //
    ALIGNMENT(
        QUALITY_CONTROL.out.samples,
        file(params.star_index),
        file(params.gtf)
    )

    //
    // SUBWORKFLOW: Analysis - RiboMetric and RiboWaltz
    // Both tools provide complementary QC and offset calculation
    //
    ANALYSIS(
        ALIGNMENT.out.transcriptome_bam,
        params.ribometric_annotation ? file(params.ribometric_annotation) : null,
        file(params.gtf)
    )

    //
    // SUBWORKFLOW: Post-processing (filter, BEDgraph, BigWig)
    //
    POST_PROCESSING(
        ALIGNMENT.out.genome_bam,
        ANALYSIS.out.offsets,
        file(params.chrom_sizes_file)
    )
}

/*
========================================================================================
    THE END
========================================================================================
*/