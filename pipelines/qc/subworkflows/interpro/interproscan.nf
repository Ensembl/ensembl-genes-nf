
include { INTERPRO_RUN } from '../../modules/interpro/run_interpro.nf'
include { PROCESS_INTERPRO } from '../../modules/interpro/interpro_stats.nf'

workflow INTERPRO_SCAN {

    take:
        protein_ch
        database
        data_file_path
        ensembl_genes_repo
        interpro_parser

    main:
        interpro_run = INTERPRO_RUN(
            protein_ch,
            database,
            data_file_path
        )

        interpro_stats_input = interpro_run.stats_txt.join(protein_ch)

        interpro_stats = PROCESS_INTERPRO(
            interpro_stats_input,
            ensembl_genes_repo,
            interpro_parser
        )

    emit:
        stats_tsv = interpro_stats.interpro_stats_tsv
        versions = interpro_run.versions.mix(interpro_stats.versions)
}
