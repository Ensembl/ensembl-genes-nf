nextflow.enable.dsl = 2

include { sncRNA } from './subworkflows/local/sncRNA.nf'

log.info """\
    s n c R N A    N F    P I P E L I N E
    =========================================
    species:     ${params.species}
    accession:   ${params.accession}
    genome file: ${params.genome_file}
    =========================================

"""



workflow {

    sncRNA(params.genome_file, params.trnascan_files, params.rfam_files)
}

