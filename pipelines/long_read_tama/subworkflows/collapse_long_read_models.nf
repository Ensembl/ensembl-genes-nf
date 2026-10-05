nextflow.enable.dsl = 2

include { SPLIT_BAM_BY_CONTIG } from '../modules/split_contigs.nf'
include { INSPECT_BAM_WORKLOAD } from '../modules/inspect_bam_workload.nf'
include { TAMA_COLLAPSE } from '../modules/tama_collapse.nf'
include { STRINGTIE2_COLLAPSE } from '../modules/stringtie2_collapse.nf'
include { STRINGTIE3_COLLAPSE } from '../modules/stringtie3_collapse.nf'
include { BAM_TO_ALIGNMENT_GTF } from '../modules/bam_to_alignment_gtf.nf'
include { TMERGE_COLLAPSE } from '../modules/tmerge_collapse.nf'
include { VALIDATE_BACKEND_BED } from '../modules/validate_backend_bed.nf'
include { AUDIT_BACKEND_SHARD_PLAN } from '../modules/audit_backend_shard_plan.nf'
include { VALIDATE_BACKEND_SHARD_PLAN } from '../modules/validate_backend_shard_plan.nf'

def backend_skip_keys(backend, explicit_skips, legacy_skips) {
    def keys = explicit_skips.findAll { value -> value.startsWith("${backend}:") }
        .collect { value -> value.substring(backend.size() + 1) }.toSet()
    if (backend == 'tama') {
        keys.addAll(legacy_skips)
    }
    return keys
}

workflow COLLAPSE_LONG_READ_MODELS {
    take:
    bam
    reference

    main:
    INSPECT_BAM_WORKLOAD(bam)
    inspected_bam = INSPECT_BAM_WORKLOAD.out.workload
    requested_backends = params.model_backend == 'all' ? ['tama', 'stringtie2', 'stringtie3', 'tmerge'] : [params.model_backend]
    split_versions = channel.empty()
    shard_plans = channel.empty()
    legacy_skips = (params.skip_tama_shards ?: '').split(',').collect { value -> value.trim() }.findAll { value -> value }.toSet()
    explicit_skips = (params.skip_model_shards ?: '').split(',').collect { value -> value.trim() }.findAll { value -> value }.toSet()
    if (params.shard_mode == 'contig') {
        SPLIT_BAM_BY_CONTIG(inspected_bam)
        split_versions = SPLIT_BAM_BY_CONTIG.out.versions
        shard_plans = SPLIT_BAM_BY_CONTIG.out.shards.map { meta, _shard_dir, manifest -> tuple(meta, manifest) }
        shard_bams = SPLIT_BAM_BY_CONTIG.out.shards.flatMap { meta, shard_dir, manifest ->
            manifest.readLines().drop(1).findAll { value -> value.trim() }.collect { line ->
                def fields = line.split('\\t', -1)
                def contig = fields[0]
                def resource_class = fields[3]
                def mapped_reads = fields[2].toLong()
                def bam_file = file("${shard_dir}/${fields[4]}")
                tuple(meta, contig, resource_class, mapped_reads, bam_file, file("${bam_file}.bai"))
            }
        }
    } else {
        shard_plans = inspected_bam.map { meta, _bam_file, _bai, workload -> tuple(meta, workload) }
        shard_bams = inspected_bam.map { meta, bam_file, bai, workload ->
            def mapped_reads = workload.readLines().drop(1).findAll { value -> value.trim() }.collect { line -> line.split('\\t', -1)[2].toLong() }.sum()
            def resource_class = mapped_reads >= params.shard_contig_reads ? 'large' : 'small'
            tuple(meta, 'whole', resource_class, mapped_reads, bam_file, bai)
        }
    }

    native_models = channel.empty()
    backend_status = channel.empty()
    status_ch = channel.empty()
    backend_shard_plans = shard_plans.flatMap { meta, plan ->
        requested_backends.collect { backend ->
            tuple(meta, backend, plan, backend_skip_keys(backend, explicit_skips, legacy_skips).toList().sort().join(','))
        }
    }
    AUDIT_BACKEND_SHARD_PLAN(backend_shard_plans)
    if (requested_backends.contains('tama')) {
        tama_bams = shard_bams.filter { meta, shard, _resource_class, _mapped_reads, _bam, _bai -> !backend_skip_keys('tama', explicit_skips, legacy_skips).contains("${meta.id}:${shard}") }
        TAMA_COLLAPSE(tama_bams, reference)
        tama_beds = TAMA_COLLAPSE.out.bed.map { meta, shard, bed -> tuple(meta, 'tama', shard, bed) }
        VALIDATE_BACKEND_BED(tama_beds)
        backend_status = VALIDATE_BACKEND_BED.out.status
        status_ch = backend_status
        native_models = native_models.mix(VALIDATE_BACKEND_BED.out.bed.map { meta, _backend, shard, bed -> tuple(meta, 'tama', shard, bed, 'bed12') })
    }
    if (requested_backends.contains('stringtie2')) {
        stringtie2_bams = shard_bams.filter { meta, shard, _resource_class, _mapped_reads, _bam, _bai -> !backend_skip_keys('stringtie2', explicit_skips, legacy_skips).contains("${meta.id}:${shard}") }
        STRINGTIE2_COLLAPSE(stringtie2_bams)
        native_models = native_models.mix(STRINGTIE2_COLLAPSE.out.gtf.map { meta, shard, gtf -> tuple(meta, 'stringtie2', shard, gtf, 'gtf') })
        status_ch = status_ch.mix(STRINGTIE2_COLLAPSE.out.status)
    }
    if (requested_backends.contains('stringtie3')) {
        stringtie3_bams = shard_bams.filter { meta, shard, _resource_class, _mapped_reads, _bam, _bai -> !backend_skip_keys('stringtie3', explicit_skips, legacy_skips).contains("${meta.id}:${shard}") }
        STRINGTIE3_COLLAPSE(stringtie3_bams)
        native_models = native_models.mix(STRINGTIE3_COLLAPSE.out.gtf.map { meta, shard, gtf -> tuple(meta, 'stringtie3', shard, gtf, 'gtf') })
        status_ch = status_ch.mix(STRINGTIE3_COLLAPSE.out.status)
    }
    if (requested_backends.contains('tmerge')) {
        tmerge_bams = shard_bams.filter { meta, shard, _resource_class, _mapped_reads, _bam, _bai -> !backend_skip_keys('tmerge', explicit_skips, legacy_skips).contains("${meta.id}:${shard}") }
        BAM_TO_ALIGNMENT_GTF(tmerge_bams)
        tmerge_gtf = BAM_TO_ALIGNMENT_GTF.out.gtf.map { meta, shard, reads_gtf -> tuple(meta, shard, reads_gtf) }
        TMERGE_COLLAPSE(tmerge_gtf)
        native_models = native_models.mix(TMERGE_COLLAPSE.out.gtf.map { meta, shard, gtf -> tuple(meta, 'tmerge', shard, gtf, 'gtf') })
        status_ch = status_ch.mix(TMERGE_COLLAPSE.out.status)
    }

    collapse_reports = INSPECT_BAM_WORKLOAD.out.workload
    collapse_reports = collapse_reports.mix(AUDIT_BACKEND_SHARD_PLAN.out.plan)
    if (requested_backends.contains('tama')) { collapse_reports = collapse_reports.mix(TAMA_COLLAPSE.out.status).mix(TAMA_COLLAPSE.out.stderr) }
    if (requested_backends.contains('stringtie2')) { collapse_reports = collapse_reports.mix(STRINGTIE2_COLLAPSE.out.status) }
    if (requested_backends.contains('stringtie3')) { collapse_reports = collapse_reports.mix(STRINGTIE3_COLLAPSE.out.status) }
    if (requested_backends.contains('tmerge')) { collapse_reports = collapse_reports.mix(TMERGE_COLLAPSE.out.status) }
    success_keys = native_models
        .map { meta, backend, shard, _model, _native_format -> "${backend}:${meta.id}:${shard}" }
        .collect()
        .ifEmpty([''])
        .map { keys -> keys.sort().join(',') }
    VALIDATE_BACKEND_SHARD_PLAN(backend_shard_plans, success_keys)
    collapse_reports = collapse_reports.mix(VALIDATE_BACKEND_SHARD_PLAN.out.status)
    version_ch = INSPECT_BAM_WORKLOAD.out.versions.mix(split_versions)
    if (requested_backends.contains('tama')) { version_ch = version_ch.mix(TAMA_COLLAPSE.out.versions) }
    if (requested_backends.contains('tama')) { version_ch = version_ch.mix(VALIDATE_BACKEND_BED.out.versions) }
    if (requested_backends.contains('stringtie2')) { version_ch = version_ch.mix(STRINGTIE2_COLLAPSE.out.versions) }
    if (requested_backends.contains('stringtie3')) { version_ch = version_ch.mix(STRINGTIE3_COLLAPSE.out.versions) }
    if (requested_backends.contains('tmerge')) { version_ch = version_ch.mix(BAM_TO_ALIGNMENT_GTF.out.versions).mix(TMERGE_COLLAPSE.out.versions) }

    emit:
    native_models
    statuses = status_ch
    shard_plans = AUDIT_BACKEND_SHARD_PLAN.out.plan
    shard_statuses = VALIDATE_BACKEND_SHARD_PLAN.out.status
    collapse_reports
    versions = version_ch
}
