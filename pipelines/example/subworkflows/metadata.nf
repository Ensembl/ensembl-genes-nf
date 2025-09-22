include { SET_MISSING_METADATA } from '../modules/set_missing_metadata.nf'
include { SET_META_CORDS } from '../modules/set_meta_coords.nf'
include { FINAL_VERIFICATION_DATACHECKS } from '../modules/metadata_datachecks.nf'
include { FIX_XREFS } from '../modules/xref_analysis_fix.nf'
include { TIDY_DUPLICATE_COORD_SYSTEMS } from '../modules/tidy_duplicate_coord_systems.nf'
include { CLEAN_ORPHANED_METADATA } from '../modules/clean_orphan_feature_metadata.nf'

workflow METADATA {

    take:
    core_info
    skip_dc

    main:
    SET_MISSING_METADATA            (core_info)
    SET_META_CORDS                  (core_info)
    FIX_XREFS                       (core_info)
    TIDY_DUPLICATE_COORD_SYSTEMS    (core_info)
    CLEAN_ORPHANED_METADATA         (core_info)


    if (!skip_dc) {
        log.info("Running final verification datachecks")
        FINAL_VERIFICATION_DATACHECKS(
            core_info,
            SET_MISSING_METADATA.out.results,
            SET_META_CORDS.out.results
        )
    } else {
        log.info("Skipping final verification datachecks")
    }


}