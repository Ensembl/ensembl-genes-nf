#!/usr/bin/env nextflow
/*
========================================================================================
    HUMAN PANGENOME PROJECTION PIPELINE
========================================================================================
    Project reference gene annotation (e.g. GENCODE on GRCh38) onto target assemblies
    (e.g. the HPRC pangenome) via whole-genome alignment + feature-level coordinate
    projection, rescue, and validation.

    Aligner is a swappable front-end module (MINIMAP2_WGA); the hpp tool performs the
    projection. See subworkflows/pangenome_projection.nf.
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { validateParameters } from 'plugin/nf-schema'

// Validate parameters against schema on startup
validateParameters()

include { PANGENOME_PROJECTION } from './subworkflows/pangenome_projection.nf'

workflow {
    // Build the per-target channel from the sample sheet: [ meta, target_fasta ]
    targets = Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [ id: row.target_id ]
            [ meta, file(row.target_fasta) ]
        }

    PANGENOME_PROJECTION(
        targets,
        file(params.ref_fasta),
        file(params.ref_gff)
    )
}
