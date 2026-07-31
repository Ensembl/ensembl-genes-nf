include { FETCH_REFSEQ } from '../../modules/fetch_refseq.nf'
include { LOAD_REFSEQ } from '../../modules/load_refseq.nf'
include { GET_METADATA } from '../../modules/get_metadata.nf'
include { LOAD_METADATA } from '../../modules/load_metadata.nf'


workflow IMPORT_REFSEQ {

    take:
        // [meta, files]
        samples_ch
        // server settings
        db_config_ch

    main:

        FETCH_REFSEQ(samples_ch)

        LOAD_REFSEQ(
            FETCH_REFSEQ.out.files,
            db_config_ch
        )

        metadata_input = LOAD_REFSEQ.out.files.map { meta ->

                def speciesParts = meta.species.tokenize(' ')
                def speciesToken = "${speciesParts[0].toLowerCase()}_${speciesParts[1].toLowerCase()}"

                def accessionToken = meta.id
                    .toLowerCase()
                    .replace("_", "")
                    .replaceFirst(/\./, "v")

                def dbName = "${speciesToken}_${accessionToken}_rs_core_114_1"

                tuple(meta + [db_name: dbName])
            }

        GET_METADATA(
                metadata_input,
                db_config_ch
            )

        LOAD_METADATA(
            GET_METADATA.out.sql,
            db_config_ch
        )
    }
