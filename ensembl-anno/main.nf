nextflow.enable.dsl = 2

include { sncRNA } from './workflows/local/sncrna.nf'

log.info """\
    
    s n c R N A    N F    P I P E L I N E
    =========================================
    species:     ${params.species}
    accession:   ${params.accession}
    genome file: ${params.genome_file}
    =========================================

"""



workflow {

    sncRNA()
}

