#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FETCH_REFSEQ } from '../modules/fetch_refseq.nf'
include { LOAD_REFSEQ } from '../modules/load_refseq.nf'
include { GET_METADATA } from '../modules/get_metadata.nf'
include { LOAD_METADATA } from '../modules/load_metadata.nf'


workflow IMPORT_REFSEQ {

    take:
        // channel: val(meta)
        //   meta: [id: RefSeq assembly accession, species: binomial species name]
        samples_ch

        // value channel: tuple val(db_host), val(db_port), val(db_user), val(db_password)
        db_config_ch

    main:

        FETCH_REFSEQ(samples_ch)

        LOAD_REFSEQ(
            FETCH_REFSEQ.out.refseq,
            db_config_ch
        )

        metadata_input = LOAD_REFSEQ.out.loaded.map { meta, loaded_marker ->

                def speciesParts = meta.species.tokenize(' ')
                def speciesToken = "${speciesParts[0].toLowerCase()}_${speciesParts[1].toLowerCase()}"

                def accessionToken = meta.id
                    .toLowerCase()
                    .replace("_", "")
                    .replaceFirst(/\./, "v")

                def dbName = "${speciesToken}_${accessionToken}_rs_core_114_1"

                meta + [db_name: dbName]
            }

        GET_METADATA(
                metadata_input,
                db_config_ch
            )

        LOAD_METADATA(
            GET_METADATA.out.sql,
            db_config_ch
        )

        versions_ch = FETCH_REFSEQ.out.versions
            .mix(LOAD_REFSEQ.out.versions)
            .mix(GET_METADATA.out.versions)
            .mix(LOAD_METADATA.out.versions)

    emit:
        loaded_refseq = LOAD_REFSEQ.out.loaded
        metadata_sql = GET_METADATA.out.sql
        loaded_metadata = LOAD_METADATA.out.loaded
        versions = versions_ch
    }
