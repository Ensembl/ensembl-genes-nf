include { SAMTOOLS_VIEW_SHARD } from '../modules/models.nf'
include { FLAIR_JUNCTIONS; FLAIR_TRANSCRIPTOME; FLAIR_EMPTY; FLAIR_COMBINE } from '../modules/optional_backends.nf'

workflow RUN_FLAIR {
    take:
    shards
    reference
    flair_threads
    flair_args

    main:
    SAMTOOLS_VIEW_SHARD(shards)
    FLAIR_JUNCTIONS(SAMTOOLS_VIEW_SHARD.out.sam)
    flair_evidence = FLAIR_JUNCTIONS.out.junctions.join(FLAIR_JUNCTIONS.out.status, by: [0, 1])
    flair_inputs = flair_evidence.join(shards, by: [0, 1]).map { meta, shard_meta, junctions, junction_status, bam, bai -> tuple(meta, shard_meta, bam, bai, junctions, junction_status) }
    flair_nonempty = flair_inputs.filter { _meta, _shard_meta, _bam, _bai, _junctions, junction_status -> junction_status.text.contains('HAS_JUNCTION') }.map { meta, shard_meta, bam, bai, junctions, _junction_status -> tuple(meta, shard_meta, bam, bai, junctions) }
    flair_empty = FLAIR_JUNCTIONS.out.status.filter { _meta, _shard_meta, status -> status.text.contains('EMPTY') }.map { meta, shard_meta, status -> tuple(meta, shard_meta, status) }
    FLAIR_TRANSCRIPTOME(flair_nonempty, reference, flair_threads, flair_args)
    FLAIR_EMPTY(flair_empty)
    flair_native = FLAIR_TRANSCRIPTOME.out.native_models.mix(FLAIR_EMPTY.out.native_models)
    flair_run_combined = flair_native
        .map { meta, shard_meta, _backend, product, _status -> tuple(meta.run_id, meta, product) }
        .groupTuple(by: 0)
        .map { _run_id, metas, gtfs -> tuple(metas[0], 'flair', gtfs) }
    FLAIR_COMBINE(flair_run_combined)
    emit:
    native_models = flair_native.map { meta, shard_meta, backend, product, status -> tuple(meta, [backend: backend, scope: 'run', run_id_or_cohort_id: meta.run_id, shard_id_or_merged: shard_meta.shard_id, native_format: 'gtf', model_count: 0, status: status.text.contains('EMPTY') ? 'EMPTY' : 'SUCCESS'], product, status) }
    combined = FLAIR_COMBINE.out.native_models.map { meta, backend, product, status ->
        tuple(meta, [backend: backend, scope: 'run', run_id_or_cohort_id: meta.run_id, shard_id_or_merged: 'merged', native_format: 'gtf', model_count: 0, status: status.text.contains('EMPTY') ? 'EMPTY' : 'SUCCESS'], product, status)
    }
    versions = SAMTOOLS_VIEW_SHARD.out.versions.mix(FLAIR_JUNCTIONS.out.versions).mix(FLAIR_TRANSCRIPTOME.out.versions).mix(FLAIR_EMPTY.out.versions).mix(FLAIR_COMBINE.out.versions)
}
