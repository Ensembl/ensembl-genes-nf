nextflow.enable.dsl=2

include { ENA_SUBMIT_WORKFLOW } from './workflows/ena_submit.nf'

workflow {
    ENA_SUBMIT_WORKFLOW()
}
