include { EVOLUTIONARY_HISTORY } from '../modules/evolutionary_history.nf'
include { TYPED_AXIS as EVOLUTIONARY_HISTORY_UNINFORMATIVE } from '../modules/typed_axis.nf'

workflow EVOLUTIONARY_HISTORY_WORKFLOW {
    take:
    instances // [meta(instance), reconciled instances]

    main:
    if (params.evolutionary_candidates) {
        evidence = channel.value(file(params.evolutionary_candidates, checkIfExists: true))
        EVOLUTIONARY_HISTORY(instances, evidence)
        axes = EVOLUTIONARY_HISTORY.out.axis.map { meta, file -> tuple(meta, 'evolutionary_history', file) }
        versions = EVOLUTIONARY_HISTORY.out.versions
    } else {
        requests = instances.map { meta, file -> tuple(meta, 'evolutionary_history', file) }
        EVOLUTIONARY_HISTORY_UNINFORMATIVE(requests)
        axes = EVOLUTIONARY_HISTORY_UNINFORMATIVE.out.axis
        versions = EVOLUTIONARY_HISTORY_UNINFORMATIVE.out.versions
    }

    emit:
    axes = axes
    versions = versions
}
