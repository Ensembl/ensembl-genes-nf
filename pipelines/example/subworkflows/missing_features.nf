include { FIND_MISSING_TRANSLATIONS } from '../modules/find_missing_translations'
include { DELETE_TRANSCRIPTS } from '../modules/delete_transcripts'

workflow MISSING_TRANSLATIONS {

    take:
    core_info 

    main:
    FIND_MISSING_TRANSLATIONS   (core_info)
    // DELETE_TRANSCRIPTS          (FIND_MISSING_TRANSLATIONS.out.results, $params.db_pass)

    // emit:
    // results = DELETE_TRANSCRIPTS.out.results
}