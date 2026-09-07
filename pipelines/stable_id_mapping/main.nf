#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RESOLVE_AND_ROUTE } from './subworkflows/local/resolve_and_route.nf'
include { MAPPING_BRANCH } from './subworkflows/local/mapping_branch.nf'
include { REASSIGNMENT_BRANCH } from './subworkflows/local/reassignment_branch.nf'


workflow {
    RESOLVE_AND_ROUTE()

    MAPPING_BRANCH(
        RESOLVE_AND_ROUTE.out.mapping
    )

    REASSIGNMENT_BRANCH(
        RESOLVE_AND_ROUTE.out.reassignment
    )
}
