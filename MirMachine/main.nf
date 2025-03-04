nextflow.enable.dsl = 2

include { mirMachine } from './subworkflows/local/mirMachine.nf'

log.info """\

    m i R N A    N F    P I P E L I N E
    =========================================
    input: ${params.input}
    =========================================
"""

workflow {
    Channel
        .fromPath(params.input)
        .splitCsv(header:true, sep:'\t')
        .map { row -> 
            def scientific_name = row.Scientific_name ?: row.'Scientific name'
            def accession = row.Accession
            def meta = [id: accession, species: scientific_name]
            tuple(meta, scientific_name, accession)
        }
        .set { input_ch }

    mirMachine(input_ch, params.fasta_dir)
}

workflow.onComplete {
    log.info "Pipeline completed at: ${new Date().format('dd-MM-yyyy HH:mm:ss')}"
}
