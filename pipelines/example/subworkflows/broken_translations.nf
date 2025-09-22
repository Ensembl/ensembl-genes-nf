

include { FIND_BROKEN_TRANSLATIONS } from '../modules/find_broken_translations'
include { DELETE_TRANSCRIPTS } from '../modules/delete_transcripts'
include { ENSURE_CANONICAL_TX } from '../modules/ensure_canonical_tx'

workflow BROKEN_TRANSLATIONS {

    take:
    core_info
    db_pass

    main:
    FIND_BROKEN_TRANSLATIONS    (core_info)
    DELETE_TRANSCRIPTS          (FIND_BROKEN_TRANSLATIONS.out.results, db_pass)
    ENSURE_CANONICAL_TX         (core_info, DELETE_TRANSCRIPTS.out.logs)

    // emit:
    // results = DELETE_TRANSCRIPTS.out.results
}