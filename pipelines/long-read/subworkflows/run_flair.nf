include { SAMTOOLS_VIEW_SHARD } from '../modules/models.nf'
include { FLAIR_JUNCTIONS; FLAIR_TRANSCRIPTOME; FINALIZE_FLAIR_GTF; COLLECT_FLAIR_PRODUCTS; WRITE_FLAIR_MANIFEST; FLAIR_EMPTY; FLAIR_COMBINE_GTF; FLAIR_COMBINE_PRODUCTS; WRITE_FLAIR_COMBINED_MANIFEST } from '../modules/optional_backends.nf'

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
    FINALIZE_FLAIR_GTF(FLAIR_TRANSCRIPTOME.out.raw)
    COLLECT_FLAIR_PRODUCTS(FLAIR_TRANSCRIPTOME.out.raw)
    flair_manifest_inputs = FINALIZE_FLAIR_GTF.out.native_models
        .join(COLLECT_FLAIR_PRODUCTS.out.products, by: [0, 1])
        .join(flair_nonempty.map { meta, shard_meta, bam, _bai, junctions -> tuple(meta, shard_meta, bam, junctions) }, by: [0, 1])
        .map { meta, shard_meta, native_gtf, status, native_products, bam, junctions -> tuple(meta, shard_meta, native_products, status, bam, junctions) }
    WRITE_FLAIR_MANIFEST(flair_manifest_inputs, reference, flair_args)
    FLAIR_EMPTY(flair_empty)
    flair_native = FINALIZE_FLAIR_GTF.out.native_models.mix(FLAIR_EMPTY.out.native_models)
    flair_run_combined = flair_native
        .map { meta, shard_meta, _backend, product, _status -> tuple(meta.run_id, meta, product) }
        .groupTuple(by: 0)
        .map { _run_id, metas, gtfs -> tuple(metas[0], 'flair', gtfs) }
    FLAIR_COMBINE_GTF(flair_run_combined)
    FLAIR_COMBINE_PRODUCTS(flair_run_combined)
    combined_manifest_inputs = FLAIR_COMBINE_GTF.out.native_models
        .join(FLAIR_COMBINE_PRODUCTS.out.products, by: [0])
        .map { meta, _backend, _native_gtf, status, gtfs, native_products -> tuple(meta, native_products, status, gtfs) }
    WRITE_FLAIR_COMBINED_MANIFEST(combined_manifest_inputs)
    emit:
    native_models = flair_native.map { meta, shard_meta, backend, product, status -> tuple(meta, [backend: backend, scope: 'run', run_id_or_cohort_id: meta.run_id, shard_id_or_merged: shard_meta.shard_id, native_format: 'gtf', model_count: 0, status: status.text.contains('EMPTY') ? 'EMPTY' : 'SUCCESS'], product, status) }
    combined = FLAIR_COMBINE_GTF.out.native_models.map { meta, backend, product, status ->
        tuple(meta, [backend: backend, scope: 'run', run_id_or_cohort_id: meta.run_id, shard_id_or_merged: 'merged', native_format: 'gtf', model_count: 0, status: status.text.contains('EMPTY') ? 'EMPTY' : 'SUCCESS'], product, status)
    }
    versions = SAMTOOLS_VIEW_SHARD.out.versions.mix(FLAIR_JUNCTIONS.out.versions).mix(FLAIR_TRANSCRIPTOME.out.versions).mix(FINALIZE_FLAIR_GTF.out.versions).mix(COLLECT_FLAIR_PRODUCTS.out.versions).mix(WRITE_FLAIR_MANIFEST.out.versions).mix(FLAIR_EMPTY.out.versions).mix(FLAIR_COMBINE_GTF.out.versions).mix(FLAIR_COMBINE_PRODUCTS.out.versions).mix(WRITE_FLAIR_COMBINED_MANIFEST.out.versions)
}
