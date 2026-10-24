nextflow.enable.dsl = 2

include { FLAIR_JUNCTIONS; FLAIR_TRANSCRIPTOME; FLAIR_COMBINE; FLAIR_COMBINE as FLAIR_COMBINE_COHORT } from '../modules/flair.nf'
include { AUDIT_NATIVE_SHARD_MODELS } from '../modules/audit_native_shard_models.nf'
include { GTF_TO_BED12 } from '../modules/gtf_to_bed12.nf'

workflow RUN_FLAIR {
    take:
    flair_shards
    reference

    main:
    // FLAIR transcriptome runs independently on each existing BAM shard.
    // Non-empty shard products are recombined per accession, then across the
    // cohort, so FLAIR's native merge behavior remains in use at both levels.
    FLAIR_JUNCTIONS(flair_shards)
    flair_inputs = flair_shards
        .map { meta, shard, resource_class, mapped_reads, bam, bai ->
            tuple("${meta.id}:${shard}", meta, shard, resource_class, mapped_reads, bam, bai)
        }
        .join(FLAIR_JUNCTIONS.out.bed.map { meta, shard, junctions -> tuple("${meta.id}:${shard}", junctions) })
        .map { _key, meta, shard, resource_class, mapped_reads, bam, bai, junctions ->
            tuple(meta, shard, resource_class, mapped_reads, bam, bai, junctions)
        }
    FLAIR_TRANSCRIPTOME(flair_inputs, reference)

    audit_inputs = FLAIR_TRANSCRIPTOME.out.gtf.map { meta, shard, gtf ->
        tuple(meta, 'flair', shard, gtf, 'gtf')
    }
    AUDIT_NATIVE_SHARD_MODELS(audit_inputs)
    nonempty_shard_keys = AUDIT_NATIVE_SHARD_MODELS.out.mergeable
        .map { meta, _backend, shard, _model, _format -> "${meta.id}:${shard}" }
    shard_archives = FLAIR_TRANSCRIPTOME.out.products
        .map { meta, shard, archive -> tuple("${meta.id}:${shard}", meta, archive) }
        .join(nonempty_shard_keys.map { key -> tuple(key, true) })
        .map { _key, meta, archive, _present -> tuple(meta.id, meta, archive) }
    accession_input = shard_archives.groupTuple()
        .map { accession, metas, archives -> tuple(metas[0] + [scope: 'accession'], 'accession', archives) }
    FLAIR_COMBINE(accession_input)

    accession_products = FLAIR_COMBINE.out.products
    cohort_input = accession_products
        .map { _meta, archive -> archive }
        .collect()
        .map { archives -> tuple([id: params.cohort_id, scope: 'cohort'], 'cohort', archives) }
    FLAIR_COMBINE_COHORT(cohort_input)

    accession_gtf = FLAIR_COMBINE.out.gtf.map { meta, gtf -> tuple(meta, 'flair', 'accession', gtf) }
    cohort_gtf = FLAIR_COMBINE_COHORT.out.gtf.map { meta, gtf -> tuple(meta, 'flair', 'cohort', gtf) }
    GTF_TO_BED12(accession_gtf.mix(cohort_gtf))

    emit:
    products = accession_products.mix(FLAIR_COMBINE_COHORT.out.products)
    bed = GTF_TO_BED12.out.bed.filter { _meta, _backend, scope, _bed -> scope == 'cohort' }
        .map { meta, _backend, scope, bed -> tuple(meta, 'flair', scope, bed) }
    accession_bed = GTF_TO_BED12.out.bed.filter { _meta, _backend, scope, _bed -> scope == 'accession' }
        .map { meta, _backend, scope, bed -> tuple(meta, 'flair', scope, bed) }
    reports = accession_products.mix(FLAIR_COMBINE_COHORT.out.products)
    model_counts = AUDIT_NATIVE_SHARD_MODELS.out.counts
    versions = FLAIR_TRANSCRIPTOME.out.versions
        .mix(FLAIR_JUNCTIONS.out.versions)
        .mix(AUDIT_NATIVE_SHARD_MODELS.out.versions)
        .mix(FLAIR_COMBINE.out.versions)
        .mix(FLAIR_COMBINE_COHORT.out.versions)
        .mix(GTF_TO_BED12.out.versions)
}
