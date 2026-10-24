nextflow.enable.dsl = 2

include { SPLIT_BAM_BY_CONTIG } from '../modules/split_contigs.nf'
include { INSPECT_BAM_WORKLOAD } from '../modules/inspect_bam_workload.nf'

workflow PREPARE_CONTIG_SHARDS {
    take:
    bam

    main:
    if (params.shard_mode == 'contig') {
        INSPECT_BAM_WORKLOAD(bam)
        SPLIT_BAM_BY_CONTIG(INSPECT_BAM_WORKLOAD.out.workload)
        shards = SPLIT_BAM_BY_CONTIG.out.shards.flatMap { meta, shard_dir, manifest ->
            manifest.readLines().drop(1).findAll { value -> value.trim() }.collect { line ->
                def fields = line.split('\\t', -1)
                def shard = fields[0]
                def resource_class = fields[3]
                def mapped_reads = fields[2].toLong()
                def bam_file = file("${shard_dir}/${fields[4]}")
                tuple(meta, shard, resource_class, mapped_reads, bam_file, file("${bam_file}.bai"))
            }
        }
        shard_versions = INSPECT_BAM_WORKLOAD.out.versions.mix(SPLIT_BAM_BY_CONTIG.out.versions)
    } else {
        INSPECT_BAM_WORKLOAD(bam)
        shards = INSPECT_BAM_WORKLOAD.out.workload.map { meta, bam_file, bai, workload ->
            def mapped_reads = workload.readLines().drop(1).findAll { value -> value.trim() }
                .collect { line -> line.split('\\t', -1)[2].toLong() }.sum()
            def resource_class = mapped_reads >= params.shard_contig_reads ? 'large' : 'small'
            tuple(meta, 'whole', resource_class, mapped_reads, bam_file, bai)
        }
        shard_versions = INSPECT_BAM_WORKLOAD.out.versions
    }

    emit:
    shards
    versions = shard_versions
}
