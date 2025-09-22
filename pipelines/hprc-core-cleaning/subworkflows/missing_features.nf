include { FIND_MISSING_TRANSLATIONS } from '../modules/find_missing_translations'
include { DELETE_TRANSCRIPTS } from '../modules/delete_transcripts'
include { ENSURE_CANONICAL_TX } from '../modules/ensure_canonical_tx'

workflow MISSING_TRANSLATIONS {

    take:
    core_info 
    db_pass

    main:
    FIND_MISSING_TRANSLATIONS   (core_info)
    DELETE_TRANSCRIPTS          (FIND_MISSING_TRANSLATIONS.out.results, db_pass)
    ENSURE_CANONICAL_TX         (core_info, DELETE_TRANSCRIPTS.out.logs)

    // emit:
    // results = DELETE_TRANSCRIPTS.out.results
}