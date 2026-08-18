include { FETCH_REFSEQ } from '../modules/fetch_refseq.nf'
include { LOAD_REFSEQ } from '../modules/load_refseq.nf'

workflow IMPORT_REFSEQ_TO_CORE {

    take:
        // channel: val(meta)
        //   meta: [id: RefSeq assembly accession, species: binomial species name]
        samples_ch

        // value channel: tuple val(db_host), val(db_port), val(db_user),
        // val(db_password), val(db_read_user)
        db_config_ch

    main:

        db_write_config_ch = db_config_ch.map { db_host, db_port, db_user, db_password, db_read_user ->
            tuple(db_host, db_port, db_user, db_password)
        }

        FETCH_REFSEQ(samples_ch)

        LOAD_REFSEQ(
            FETCH_REFSEQ.out.refseq,
            db_write_config_ch
        )

        metadata_input = LOAD_REFSEQ.out.loaded.map { meta, loaded_marker ->

                def speciesParts = meta.species.tokenize(' ')
                def speciesToken = "${speciesParts[0].toLowerCase()}_${speciesParts[1].toLowerCase()}"

                def accessionToken = meta.id
                    .toLowerCase()
                    .replace("_", "")
                    .replaceFirst(/\./, "v")

                def dbName = "${speciesToken}_${accessionToken}_rs_core_114_1"

                tuple(meta + [db_name: dbName], loaded_marker)
            }

        versions_ch = FETCH_REFSEQ.out.versions
            .mix(LOAD_REFSEQ.out.versions)


    emit:
        loaded_refseq = LOAD_REFSEQ.out.loaded
        metadata_input
        versions = versions_ch
    }
