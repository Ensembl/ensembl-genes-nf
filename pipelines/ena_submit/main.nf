nextflow.enable.dsl=2

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include { ENA_SUBMIT_WORKFLOW } from './workflows/ena_submit.nf'

workflow {
    validateParameters()
    log.info(paramsSummaryLog(workflow))
    ENA_SUBMIT_WORKFLOW()
}
