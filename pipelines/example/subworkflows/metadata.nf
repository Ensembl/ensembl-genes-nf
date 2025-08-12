include { SET_MISSING_METADATA } from '../modules/set_missing_metadata.nf'
include { FINAL_VERIFICATION_DATACHECKS } from '../modules/metadata_datachecks.nf'

workflow METADATA {

    take:
    core_info 

    main:
    SET_MISSING_METADATA            (core_info)
    FINAL_VERIFICATION_DATACHECKS   (core_info)

}