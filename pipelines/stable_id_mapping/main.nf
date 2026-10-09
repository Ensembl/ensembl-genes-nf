#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RESOLVE_AND_ROUTE } from './subworkflows/resolve_and_route.nf'
include { MAPPING_BRANCH } from './subworkflows/mapping_branch.nf'
include { REASSIGNMENT_BRANCH } from './subworkflows/reassignment_branch.nf'


workflow {
    RESOLVE_AND_ROUTE()

    MAPPING_BRANCH(
        RESOLVE_AND_ROUTE.out.mapping,
        RESOLVE_AND_ROUTE.out.mapping_metadata
    )

    REASSIGNMENT_BRANCH(
        RESOLVE_AND_ROUTE.out.reassignment
    )
}
