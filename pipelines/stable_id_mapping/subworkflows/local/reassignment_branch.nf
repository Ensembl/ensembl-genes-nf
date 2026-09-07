include { REASSIGN_STABLE_IDS } from '../../modules/local/reassign_stable_ids.nf'


workflow REASSIGNMENT_BRANCH {
    take:
    reassignment_ch

    main:
    REASSIGN_STABLE_IDS(reassignment_ch)

    emit:
    reassignment = REASSIGN_STABLE_IDS.out.reassignment
}
