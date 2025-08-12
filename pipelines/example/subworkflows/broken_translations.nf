

include { FIND_BROKEN_TRANSLATIONS } from '../modules/find_broken_translations'
include { DELETE_TRANSCRIPTS } from '../modules/delete_transcripts'

workflow BROKEN_TRANSLATIONS {

    take:
    core_info 

    main:
    FIND_BROKEN_TRANSLATIONS    (core_info)
    // DELETE_TRANSCRIPTS          (FIND_BROKEN_TRANSLATIONS.out.results, $params.db_pass)

    // emit:
    // results = DELETE_TRANSCRIPTS.out.results
}