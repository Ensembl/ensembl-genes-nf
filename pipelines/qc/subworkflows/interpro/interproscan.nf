
include { INTERPRO_RUN } from '../../modules/interpro/run_interpro.nf'
include { PROCESS_INTERPRO } from '../../modules/interpro/interpro_stats.nf'

workflow INTERPRO_SCAN {

    take:
        // [meta, protein]
        protein_ch
        //optional path to interpro database file
        data_file_path
        // optional database to run against - defaults to pfam
        database
        // path to ensembl-genes checkout
        ensembl_genes_repo

    main:
        interpro_tsv = INTERPRO_RUN(
            protein_ch,
            data_file_path,
            database
        )

        stats_tsv = PROCESS_INTERPRO(
            interpro_tsv,
            protein_ch
        )

    emit:
        stats_tsv
}
