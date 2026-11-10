include { PREPARE_BACKEND_MERGE; MERGE_DIRECT_GTF; STRINGTIE2_MERGE; STRINGTIE3_MERGE; TMERGE_MERGE; EMIT_EMPTY_MERGED_GTF; WRITE_MERGE_FILELIST } from '../modules/merge.nf'

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

    merge_inputs = run_grouped
        .map { meta, backend, products -> tuple(meta, backend, 'run', products, backend == 'tmerge' ? 'coordinates' : 'name') }
        .mix(cohort_grouped.map { meta, backend, products -> tuple(meta, backend, 'cohort', products, backend == 'tmerge' ? 'coordinates' : 'name') })
    PREPARE_BACKEND_MERGE(merge_inputs)
    prepared = PREPARE_BACKEND_MERGE.out.prepared

    direct_inputs = prepared
        .filter { _meta, backend, _scope, _gtfs, _merge_input, status, _merge_status, _skipped -> !(backend in ['stringtie2', 'stringtie3', 'tmerge']) && status.text.contains('COMPLETE') }
        .map { meta, backend, scope, _gtfs, merge_input, status, merge_status, skipped -> tuple(meta, backend, scope, merge_input, status, merge_status, skipped) }
    MERGE_DIRECT_GTF(direct_inputs)

    stringtie2_inputs = prepared
        .filter { _meta, backend, _scope, _gtfs, _merge_input, status, _merge_status, _skipped -> backend == 'stringtie2' && status.text.contains('COMPLETE') }
        .map { meta, backend, scope, gtfs, merge_input, status, _merge_status, _skipped -> tuple(meta, backend, scope, gtfs, merge_input, status) }
    STRINGTIE2_MERGE(stringtie2_inputs, stringtie2_merge_args)

    stringtie3_inputs = prepared
        .filter { _meta, backend, _scope, _gtfs, _merge_input, status, _merge_status, _skipped -> backend == 'stringtie3' && status.text.contains('COMPLETE') }
        .map { meta, backend, scope, gtfs, merge_input, status, _merge_status, _skipped -> tuple(meta, backend, scope, gtfs, merge_input, status) }
    STRINGTIE3_MERGE(stringtie3_inputs, stringtie3_merge_args)

    tmerge_inputs = prepared
        .filter { _meta, backend, _scope, _gtfs, _merge_input, status, _merge_status, _skipped -> backend == 'tmerge' && status.text.contains('COMPLETE') }
        .map { meta, backend, scope, _gtfs, merge_input, status, _merge_status, _skipped -> tuple(meta, backend, scope, merge_input, status) }
    TMERGE_MERGE(tmerge_inputs, tmerge_args)

    empty_inputs = prepared
        .filter { _meta, _backend, _scope, _gtfs, _merge_input, status, _merge_status, _skipped -> status.text.contains('EMPTY') }
        .map { meta, backend, scope, _gtfs, merge_input, status, merge_status, skipped -> tuple(meta, backend, scope, merge_input, status, merge_status, skipped) }
    EMIT_EMPTY_MERGED_GTF(empty_inputs)

    WRITE_MERGE_FILELIST(prepared.map { meta, backend, scope, gtfs, _merge_input, status, _merge_status, _skipped -> tuple(meta, backend, scope, gtfs, status) })

    merged_gtfs = MERGE_DIRECT_GTF.out.merged
        .mix(STRINGTIE2_MERGE.out.merged)
        .mix(STRINGTIE3_MERGE.out.merged)
        .mix(TMERGE_MERGE.out.merged)
        .mix(EMIT_EMPTY_MERGED_GTF.out.merged)
    merged = merged_gtfs
        .join(WRITE_MERGE_FILELIST.out.filelist, by: [0, 1, 2])
        .map { meta, backend, scope, merged_gtf, filelist, checksum -> tuple(meta, backend, scope, merged_gtf, filelist, checksum) }

    emit:
    merged
    versions = STRINGTIE2_MERGE.out.versions
        .mix(STRINGTIE3_MERGE.out.versions)
        .mix(TMERGE_MERGE.out.versions)
}
