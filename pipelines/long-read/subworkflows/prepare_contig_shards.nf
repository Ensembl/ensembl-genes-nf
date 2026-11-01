include { SAMTOOLS_IDXSTATS } from '../modules/alignment.nf'
include { INSPECT_BAM_WORKLOAD; PLAN_CONTIG_SHARDS; SPLIT_BAM_BY_CONTIG; INDEX_CONTIG_SHARD; VALIDATE_CONTIG_SHARDS; COUNT_SOURCE_COVERAGE; COUNT_CONTIG_SHARD_COVERAGE; COLLECT_SHARD_COVERAGE_COUNTS; VALIDATE_SHARD_COVERAGE } from '../modules/sharding.nf'

workflow PREPARE_CONTIG_SHARDS {
    take:
    alignment
    shard_mode
    shard_contig_reads
    shard_very_large_reads

    main:
    SAMTOOLS_IDXSTATS(alignment)
    INSPECT_BAM_WORKLOAD(SAMTOOLS_IDXSTATS.out.idxstats)
    PLAN_CONTIG_SHARDS(INSPECT_BAM_WORKLOAD.out.workload, shard_mode, shard_contig_reads, shard_very_large_reads)
    plan_rows = PLAN_CONTIG_SHARDS.out.plan.flatMap { meta, bam, bai, plan, _summary, _status ->
        plan.readLines().drop(1).findAll { line -> line.trim() && line.split('\\t', -1)[4] == 'true' }.collect { line ->
            def fields = line.split('\\t', -1)
            tuple(meta, [shard_id: fields[0], contig: fields[1], reference_length: fields[2] ? fields[2].toLong() : 0, mapped_reads: fields[3].toLong(), resource_class: fields[5], source_run_id: fields[6], shard_status: fields[7]], bam, bai)
        }
    }
    SPLIT_BAM_BY_CONTIG(plan_rows)
    INDEX_CONTIG_SHARD(SPLIT_BAM_BY_CONTIG.out.shards)
    VALIDATE_CONTIG_SHARDS(INDEX_CONTIG_SHARD.out.indexed)
    coverage_inputs = PLAN_CONTIG_SHARDS.out.plan
        .filter { _meta, _bam, _bai, _plan, summary, _status -> summary.readLines()[1].split('\\t')[1].toInteger() > 0 }
        .map { meta, bam, _bai, plan, _summary, _status -> tuple(meta.run_id, meta, bam, plan) }
        .join(VALIDATE_CONTIG_SHARDS.out.validated.map { meta, _shard_meta, bam, _bai, _validation -> tuple(meta.run_id, bam) }.groupTuple(), by: 0)
        .map { _run_id, meta, source_bam, plan, shard_bams -> tuple(meta, source_bam, plan, shard_bams) }
    COUNT_SOURCE_COVERAGE(coverage_inputs.map { meta, source_bam, plan, _shard_bams -> tuple(meta, source_bam, plan) })
    COUNT_CONTIG_SHARD_COVERAGE(VALIDATE_CONTIG_SHARDS.out.validated.map { meta, shard_meta, bam, _bai, _validation -> tuple(meta, shard_meta, bam) })
    coverage_counts = COUNT_SOURCE_COVERAGE.out.source_count
        .map { meta, plan, source_count -> tuple(meta.run_id, meta, plan, source_count) }
        .join(COUNT_CONTIG_SHARD_COVERAGE.out.shard_count.map { meta, _shard_meta, _bam, count -> tuple(meta.run_id, count) }.groupTuple(), by: 0)
        .join(VALIDATE_CONTIG_SHARDS.out.validated.map { meta, _shard_meta, bam, _bai, _validation -> tuple(meta.run_id, bam) }.groupTuple(), by: 0)
        .map { _run_id, meta, plan, source_count, shard_counts, shard_bams -> tuple(meta, plan, source_count, shard_counts, shard_bams) }
    COLLECT_SHARD_COVERAGE_COUNTS(coverage_counts)
    VALIDATE_SHARD_COVERAGE(COLLECT_SHARD_COVERAGE_COUNTS.out.counts)

    emit:
    workload = INSPECT_BAM_WORKLOAD.out.workload
    plans = PLAN_CONTIG_SHARDS.out.plan.map { meta, bam, bai, plan, _summary, _status -> tuple(meta, bam, bai, plan) }
    summaries = PLAN_CONTIG_SHARDS.out.plan.map { meta, _bam, _bai, _plan, summary, _status -> tuple(meta, summary.readLines()[1].split('\\t')[1].toInteger()) }
    shards = VALIDATE_CONTIG_SHARDS.out.validated.map { meta, shard_meta, bam, bai, validation -> tuple(meta, shard_meta, bam, bai) }
    status = SPLIT_BAM_BY_CONTIG.out.status
    versions = SAMTOOLS_IDXSTATS.out.versions.mix(INSPECT_BAM_WORKLOAD.out.versions).mix(PLAN_CONTIG_SHARDS.out.versions).mix(SPLIT_BAM_BY_CONTIG.out.versions).mix(INDEX_CONTIG_SHARD.out.versions).mix(VALIDATE_CONTIG_SHARDS.out.versions).mix(COUNT_SOURCE_COVERAGE.out.versions).mix(COUNT_CONTIG_SHARD_COVERAGE.out.versions).mix(COLLECT_SHARD_COVERAGE_COUNTS.out.versions).mix(VALIDATE_SHARD_COVERAGE.out.versions)
}
