include { TRANSLATION_JOIN } from '../modules/translation_join.nf'

workflow TRANSLATION_RECONCILIATION {
    take:
    instances // [meta, instances.jsonl]
    verdicts  // [meta, translonscorer verdict TSV/JSONL]

    main:
    TRANSLATION_JOIN(instances, verdicts)

    emit:
    instances = TRANSLATION_JOIN.out.instances
    versions = TRANSLATION_JOIN.out.versions
}
