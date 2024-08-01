nextflow.enable.dsl = 2

include { bam_processing } from "${projectDir}/subworkflows/local/bam_processing.nf"

log.info """\
    B A M s    N F    P I P E L I N E
    =========================================
    species     : ${params.species}
    assembly    : ${params.assembly}
    parent_dir  : ${params.parent_dir}
    =========================================
"""

workflow {
    bam_processing(params.species, params.assembly, params.parent_dir)
}

workflow.onComplete {
    log.info "Pipeline completed at: ${new Date().format('yyyy-MM-dd HH:mm:ss')}"
}
