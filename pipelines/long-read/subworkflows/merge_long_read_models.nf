include { TAMA_MERGE_ACCESSION as TAMA_MERGE_RUN; TAMA_MERGE as TAMA_MERGE_COHORT } from '../modules/merge.nf'
include { STRINGTIE2_MERGE as STRINGTIE2_MERGE_RUN; STRINGTIE2_MERGE as STRINGTIE2_MERGE_COHORT } from '../modules/merge.nf'
include { STRINGTIE3_MERGE as STRINGTIE3_MERGE_RUN; STRINGTIE3_MERGE as STRINGTIE3_MERGE_COHORT } from '../modules/merge.nf'
include { TMERGE_MERGE_ACCESSION as TMERGE_MERGE_RUN; TMERGE_NATIVE_MERGE as TMERGE_MERGE_COHORT } from '../modules/merge.nf'
include { MERGE_BACKEND_MODELS as MERGE_OTHER_RUN } from '../modules/merge_generic.nf'
include { MERGE_BACKEND_MODELS as MERGE_OTHER_COHORT } from '../modules/merge_generic.nf'

workflow MERGE_LONG_READ_MODELS {
    take:
    native_models
    cohort_id
    scope_policy
    stringtie2_merge_args
    stringtie3_merge_args
    tmerge_args

    main:
    run_grouped = native_models
            .map { meta, backend_meta, native_product, _status_report -> tuple(meta.run_id, backend_meta.backend, native_product) }
            .groupTuple(by: [0, 1])
            .map { run_id, backend, products -> tuple([id: run_id, run_id: run_id, cohort_id: cohort_id], backend, products) }
            .filter { _meta, backend, _products -> scope_policy[backend].contains('run') }
    cohort_grouped = native_models
            .map { _meta, backend_meta, native_product, _status_report -> tuple(backend_meta.backend, native_product) }
            .groupTuple(by: 0)
            .map { backend, products -> tuple([id: cohort_id, cohort_id: cohort_id], backend, products) }
            .filter { _meta, backend, _products -> scope_policy[backend].contains('cohort') }

    run_others = run_grouped.filter { _meta, backend, _products -> !(backend in ['tama', 'stringtie2', 'stringtie3', 'tmerge']) }
        .map { meta, backend, products -> tuple(meta, backend, products, 'run') }
    cohort_others = cohort_grouped.filter { _meta, backend, _products -> !(backend in ['tama', 'stringtie2', 'stringtie3', 'tmerge']) }
        .map { meta, backend, products -> tuple(meta, backend, products, 'cohort') }

    run_tama = run_grouped.filter { _meta, backend, _products -> backend == 'tama' }
    cohort_tama = cohort_grouped.filter { _meta, backend, _products -> backend == 'tama' }
    run_stringtie2 = run_grouped.filter { _meta, backend, _products -> backend == 'stringtie2' }
        .map { meta, backend, products -> tuple(meta, backend, 'run', products) }
    cohort_stringtie2 = cohort_grouped.filter { _meta, backend, _products -> backend == 'stringtie2' }
        .map { meta, backend, products -> tuple(meta, backend, 'cohort', products) }
    run_stringtie3 = run_grouped.filter { _meta, backend, _products -> backend == 'stringtie3' }
        .map { meta, backend, products -> tuple(meta, backend, 'run', products) }
    cohort_stringtie3 = cohort_grouped.filter { _meta, backend, _products -> backend == 'stringtie3' }
        .map { meta, backend, products -> tuple(meta, backend, 'cohort', products) }
    run_tmerge = run_grouped.filter { _meta, backend, _products -> backend == 'tmerge' }
    cohort_tmerge = cohort_grouped.filter { _meta, backend, _products -> backend == 'tmerge' }

    merged = channel.empty()
    versions = channel.empty()
    if (scope_policy.values().any { it.contains('run') }) {
        TAMA_MERGE_RUN(run_tama)
        STRINGTIE2_MERGE_RUN(run_stringtie2, stringtie2_merge_args)
        STRINGTIE3_MERGE_RUN(run_stringtie3, stringtie3_merge_args)
        TMERGE_MERGE_RUN(run_tmerge, tmerge_args)
        MERGE_OTHER_RUN(run_others)
        merged = merged.mix(TAMA_MERGE_RUN.out.merged)
            .mix(STRINGTIE2_MERGE_RUN.out.merged)
            .mix(STRINGTIE3_MERGE_RUN.out.merged)
            .mix(TMERGE_MERGE_RUN.out.merged)
            .mix(MERGE_OTHER_RUN.out.merged)
        versions = versions.mix(TAMA_MERGE_RUN.out.versions)
            .mix(STRINGTIE2_MERGE_RUN.out.versions)
            .mix(STRINGTIE3_MERGE_RUN.out.versions)
            .mix(TMERGE_MERGE_RUN.out.versions)
            .mix(MERGE_OTHER_RUN.out.versions)
    }
    if (scope_policy.values().any { it.contains('cohort') }) {
        TAMA_MERGE_COHORT(cohort_tama)
        STRINGTIE2_MERGE_COHORT(cohort_stringtie2, stringtie2_merge_args)
        STRINGTIE3_MERGE_COHORT(cohort_stringtie3, stringtie3_merge_args)
        TMERGE_MERGE_COHORT(cohort_tmerge, tmerge_args)
        MERGE_OTHER_COHORT(cohort_others)
        merged = merged.mix(TAMA_MERGE_COHORT.out.merged)
            .mix(STRINGTIE2_MERGE_COHORT.out.merged)
            .mix(STRINGTIE3_MERGE_COHORT.out.merged)
            .mix(TMERGE_MERGE_COHORT.out.merged)
            .mix(MERGE_OTHER_COHORT.out.merged)
        versions = versions.mix(TAMA_MERGE_COHORT.out.versions)
            .mix(STRINGTIE2_MERGE_COHORT.out.versions)
            .mix(STRINGTIE3_MERGE_COHORT.out.versions)
            .mix(TMERGE_MERGE_COHORT.out.versions)
            .mix(MERGE_OTHER_COHORT.out.versions)
    }

    emit:
    merged
    versions
}
