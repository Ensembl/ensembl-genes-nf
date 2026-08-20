include { IMPORT_REFSEQ_TO_CORE } from '../subworkflows/import_refseq_to_core.nf'
include { GET_METADATA_CORE } from '../subworkflows/get_metadata_core.nf'
include { PREPARE_STATS_INPUT } from '../modules/prepare_stats_input.nf'
include { COMBINE_STATS_INPUT } from '../modules/combine_stats_input.nf'

workflow IMPORT_REFSEQ {

    take:
        // channel: val(meta)
        //   meta: [id: RefSeq assembly accession, species: binomial species name]
        samples_ch

        // value channel: tuple val(db_host), val(db_port), val(db_user),
        // val(db_password), val(db_read_user)
        db_config_ch

    main:

        IMPORT_REFSEQ_TO_CORE(samples_ch, db_config_ch)

        GET_METADATA_CORE(
            IMPORT_REFSEQ_TO_CORE.out.metadata_input,
            db_config_ch
        )

        PREPARE_STATS_INPUT(
            IMPORT_REFSEQ_TO_CORE.out.genome
        )
        combined_stats_input = COMBINE_STATS_INPUT(PREPARE_STATS_INPUT.out.csv.collect())

        versions_ch = IMPORT_REFSEQ_TO_CORE.out.versions
            .mix(GET_METADATA_CORE.out.versions)

    emit:
        loaded_refseq = IMPORT_REFSEQ_TO_CORE.out.loaded_refseq
        stats_input = combined_stats_input.csv
        //metadata_sql = GET_METADATA_CORE.out.metadata_sql
        //loaded_metadata = GET_METADATA_CORE.out.loaded_metadata
        taxonomy_loaded = GET_METADATA_CORE.out.taxonomy_loaded
        versions = versions_ch
    }
