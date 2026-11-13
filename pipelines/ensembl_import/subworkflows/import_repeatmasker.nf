include { DOWNLOAD_REPEATMASKER } from '../modules/download_repeatmasker.nf'
include { CONVERT_REPEATMASKER } from '../modules/convert_repeatmasker.nf'
include { LOAD_REPEATMASKER } from '../modules/load_repeatmasker.nf'

workflow IMPORT_REPEATMASKER {

    take:
        // channel: tuple val(meta), path(refseq_loaded)
        loaded_refseq

        // channel: tuple val(meta), path(assembly_report)
        assembly_report

        // value channel: tuple val(db_host), val(db_port), val(db_user),
        // val(db_password), val(db_read_user)
        db_config_ch

    main:
        DOWNLOAD_REPEATMASKER(loaded_refseq)

        repeatmasker_downloaded = DOWNLOAD_REPEATMASKER.out.repeatmasker
            .map { meta, repeatmasker, repeatmasker_status ->
                tuple(meta.id, meta, repeatmasker, repeatmasker_status.text.trim() == 'available')
            }

        repeatmasker_with_report = repeatmasker_downloaded
            .filter { _key, _meta, _repeatmasker, has_repeatmasker -> has_repeatmasker }
            .map { _key, meta, repeatmasker, _has_repeatmasker -> tuple(meta, repeatmasker) }
            .join(assembly_report)

        CONVERT_REPEATMASKER(repeatmasker_with_report)

        db_write_config_ch = db_config_ch.map { db_host, db_port, db_user, db_password, _db_read_user ->
            tuple(db_host, db_port, db_user, db_password)
        }

        LOAD_REPEATMASKER(
            CONVERT_REPEATMASKER.out.gtf,
            db_write_config_ch
        )

        repeatmasker_completed = LOAD_REPEATMASKER.out.loaded
            .map { meta, _loaded_marker -> tuple(meta, 'loaded') }
            .mix(
                repeatmasker_downloaded
                    .filter { _key, _meta, _repeatmasker, has_repeatmasker -> !has_repeatmasker }
                    .map { _key, meta, _repeatmasker, _has_repeatmasker -> tuple(meta, 'missing') }
            )

        versions_ch = DOWNLOAD_REPEATMASKER.out.versions
            .mix(CONVERT_REPEATMASKER.out.versions)
            .mix(LOAD_REPEATMASKER.out.versions)

    emit:
        loaded = LOAD_REPEATMASKER.out.loaded
        completed = repeatmasker_completed
        versions = versions_ch
}
