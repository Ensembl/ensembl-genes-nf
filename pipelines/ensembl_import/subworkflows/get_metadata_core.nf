include { GET_SAMPLE_GENE } from '../modules/get_sample_gene.nf'
include { ADD_STATIC_METAKEYS } from '../modules/add_static_metakeys.nf'
include { GET_METADATA } from '../modules/get_metadata.nf'
include { LOAD_METADATA } from '../modules/load_metadata.nf'
include { LOAD_TAXONOMY } from '../modules/load_taxonomy.nf'

workflow GET_METADATA_CORE {

    take:
        // channel: val(meta)
        //   meta: [id: RefSeq assembly accession, species: binomial species name]
        metadata_input

        // value channel: tuple val(db_host), val(db_port), val(db_user),
        // val(db_password), val(db_read_user)
        db_config_ch

    main:

        db_write_config_ch = db_config_ch.map { db_host, db_port, db_user, db_password, db_read_user ->
            tuple(db_host, db_port, db_user, db_password)
        }

        GET_SAMPLE_GENE(
            metadata_input,
            db_write_config_ch
        )

        sample_gene_input = GET_SAMPLE_GENE.out.sample_gene.map { meta, sample_gene_marker ->
            meta
        }

        ADD_STATIC_METAKEYS(
            sample_gene_input,
            db_write_config_ch
        )

        metadata_post_static = ADD_STATIC_METAKEYS.out.loaded.map { meta, static_marker ->
            tuple(meta, static_marker)
        }

        //GET_METADATA(
        //    metadata_post_static,
        //    db_write_config_ch
        //)

        //LOAD_METADATA(
        //    GET_METADATA.out.sql,
        //    db_write_config_ch
        //)

        LOAD_TAXONOMY(
            metadata_post_static,
            db_config_ch
        )

        versions_ch = GET_SAMPLE_GENE.out.versions
            .mix(ADD_STATIC_METAKEYS.out.versions)
            //.mix(GET_METADATA.out.versions)
            //.mix(LOAD_METADATA.out.versions)
            .mix(LOAD_TAXONOMY.out.versions)

    emit:
        //metadata_sql = GET_METADATA.out.sql
        //loaded_metadata = LOAD_METADATA.out.loaded
        taxonomy_loaded = LOAD_TAXONOMY.out.taxonomy
        versions = versions_ch
    }
