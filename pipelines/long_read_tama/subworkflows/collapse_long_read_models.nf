nextflow.enable.dsl = 2

include { SPLIT_BAM_BY_CONTIG } from '../modules/split_contigs.nf'
include { INSPECT_BAM_WORKLOAD } from '../modules/inspect_bam_workload.nf'
include { TAMA_COLLAPSE } from '../modules/tama_collapse.nf'
include { STRINGTIE2_COLLAPSE } from '../modules/stringtie2_collapse.nf'
include { STRINGTIE3_COLLAPSE } from '../modules/stringtie3_collapse.nf'
include { BAM_TO_ALIGNMENT_GTF } from '../modules/bam_to_alignment_gtf.nf'
include { TMERGE_COLLAPSE } from '../modules/tmerge_collapse.nf'
include { VALIDATE_BACKEND_BED } from '../modules/validate_backend_bed.nf'

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
    legacy_skips = (params.skip_tama_shards ?: '').split(',').collect { it.trim() }.findAll { it }.toSet()
    explicit_skips = (params.skip_model_shards ?: '').split(',').collect { it.trim() }.findAll { it }.toSet()
    if (params.shard_mode == 'contig') {
        SPLIT_BAM_BY_CONTIG(inspected_bam)
        shard_bams = SPLIT_BAM_BY_CONTIG.out.shards.flatMap { meta, shard_dir, manifest ->
            manifest.readLines().drop(1).findAll { it.trim() }.collect { line ->
                def fields = line.split('\\t', -1)
                def contig = fields[0]
                def resource_class = fields[3]
                def mapped_reads = fields[2].toLong()
                def bam_file = file("${shard_dir}/${fields[4]}")
                tuple(meta, contig, resource_class, mapped_reads, bam_file, file("${bam_file}.bai"))
            }
        }
    } else {
        shard_bams = inspected_bam.map { meta, bam_file, bai, workload ->
            def mapped_reads = workload.readLines().drop(1).findAll { it.trim() }.collect { it.split('\\t', -1)[2].toLong() }.sum()
            def resource_class = mapped_reads >= params.shard_contig_reads ? 'large' : 'small'
            tuple(meta, 'whole', resource_class, mapped_reads, bam_file, bai)
        }
    }

    raw_beds = channel.empty()
    if (requested_backends.contains('tama')) {
        tama_bams = shard_bams.filter { meta, shard, _resource_class, _mapped_reads, _bam, _bai -> !backend_skip_keys('tama', explicit_skips, legacy_skips).contains("${meta.id}:${shard}") }
        TAMA_COLLAPSE(tama_bams, reference)
        raw_beds = raw_beds.mix(TAMA_COLLAPSE.out.bed.map { meta, shard, bed -> tuple('tama', meta, shard, bed) })
    }
    if (requested_backends.contains('stringtie2')) {
        stringtie2_bams = shard_bams.filter { meta, shard, _resource_class, _mapped_reads, _bam, _bai -> !backend_skip_keys('stringtie2', explicit_skips, legacy_skips).contains("${meta.id}:${shard}") }
        STRINGTIE2_COLLAPSE(stringtie2_bams)
        raw_beds = raw_beds.mix(STRINGTIE2_COLLAPSE.out.bed.map { meta, shard, bed -> tuple('stringtie2', meta, shard, bed) })
    }
    if (requested_backends.contains('stringtie3')) {
        stringtie3_bams = shard_bams.filter { meta, shard, _resource_class, _mapped_reads, _bam, _bai -> !backend_skip_keys('stringtie3', explicit_skips, legacy_skips).contains("${meta.id}:${shard}") }
        STRINGTIE3_COLLAPSE(stringtie3_bams)
        raw_beds = raw_beds.mix(STRINGTIE3_COLLAPSE.out.bed.map { meta, shard, bed -> tuple('stringtie3', meta, shard, bed) })
    }
    if (requested_backends.contains('tmerge')) {
        tmerge_bams = shard_bams.filter { meta, shard, _resource_class, _mapped_reads, _bam, _bai -> !backend_skip_keys('tmerge', explicit_skips, legacy_skips).contains("${meta.id}:${shard}") }
        BAM_TO_ALIGNMENT_GTF(tmerge_bams)
        tmerge_gtf = BAM_TO_ALIGNMENT_GTF.out.gtf.map { meta, shard, reads_gtf -> tuple(meta, shard, reads_gtf) }
        TMERGE_COLLAPSE(tmerge_gtf)
        raw_beds = raw_beds.mix(TMERGE_COLLAPSE.out.bed.map { meta, shard, bed -> tuple('tmerge', meta, shard, bed) })
    }

    VALIDATE_BACKEND_BED(raw_beds)
    collapse_reports = INSPECT_BAM_WORKLOAD.out.workload
    if (requested_backends.contains('tama')) { collapse_reports = collapse_reports.mix(TAMA_COLLAPSE.out.status).mix(TAMA_COLLAPSE.out.stderr) }
    version_ch = INSPECT_BAM_WORKLOAD.out.versions.mix(SPLIT_BAM_BY_CONTIG.out.versions).mix(VALIDATE_BACKEND_BED.out.versions)
    if (requested_backends.contains('tama')) { version_ch = version_ch.mix(TAMA_COLLAPSE.out.versions) }
    if (requested_backends.contains('stringtie2')) { version_ch = version_ch.mix(STRINGTIE2_COLLAPSE.out.versions) }
    if (requested_backends.contains('stringtie3')) { version_ch = version_ch.mix(STRINGTIE3_COLLAPSE.out.versions) }
    if (requested_backends.contains('tmerge')) { version_ch = version_ch.mix(BAM_TO_ALIGNMENT_GTF.out.versions).mix(TMERGE_COLLAPSE.out.versions) }

    emit:
    beds = VALIDATE_BACKEND_BED.out.bed
    statuses = VALIDATE_BACKEND_BED.out.status
    collapse_reports
    versions = version_ch
}
