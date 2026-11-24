#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RESOLVE_AND_ROUTE } from './subworkflows/resolve_and_route.nf'
include { MAPPING_BRANCH } from './subworkflows/mapping_branch.nf'
include { REASSIGNMENT_BRANCH } from './subworkflows/reassignment_branch.nf'
include { COMPARE_BRANCH } from './subworkflows/compare_branch.nf'


workflow {
    def requested_mode = params.mode?.toString()?.trim() ?: 'auto'

    if (requested_mode == 'compare') {
        if (params.samplesheet != null || params.db_name != null) {
            error "Mode 'compare' does not use --db_name or --samplesheet"
        }

        def required_inputs = [
            comparison_name: params.comparison_name,
            ref_fasta: params.ref_fasta,
            ref_gff: params.ref_gff,
            target_fasta: params.target_fasta,
            target_gff: params.target_gff,
        ]
        required_inputs.each { name, value ->
            if (value == null || !value.toString().trim()) {
                error "Mode 'compare' requires --${name}"
            }
        }

        COMPARE_BRANCH(
            Channel.of(
                tuple(
                    params.comparison_name.toString().trim(),
                    file(params.ref_fasta, checkIfExists: true),
                    file(params.ref_gff, checkIfExists: true),
                    file(params.target_fasta, checkIfExists: true),
                    file(params.target_gff, checkIfExists: true),
                )
            )
        )
    } else {
        RESOLVE_AND_ROUTE()

        MAPPING_BRANCH(
            RESOLVE_AND_ROUTE.out.mapping,
            RESOLVE_AND_ROUTE.out.mapping_metadata
        )

        REASSIGNMENT_BRANCH(
            RESOLVE_AND_ROUTE.out.reassignment
        )
    }
}
