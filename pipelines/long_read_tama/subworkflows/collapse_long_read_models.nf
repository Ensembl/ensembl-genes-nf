nextflow.enable.dsl = 2

include { SPLIT_BAM_BY_CONTIG } from '../modules/split_contigs.nf'
include { INSPECT_BAM_WORKLOAD } from '../modules/inspect_bam_workload.nf'
include { TAMA_COLLAPSE } from '../modules/tama_collapse.nf'
include { TAMA_MERGE as TAMA_MERGE_ACCESSION } from '../modules/tama_merge.nf'
include { VALIDATE_TAMA_OUTPUT } from '../modules/validate_tama_output.nf'

workflow COLLAPSE_LONG_READ_MODELS {
    take:
    bam
    reference
    inspector
    splitter
    validator
    merge_filelist_builder

    main:
    // bam: tuple val(meta), path(sorted.bam), path(sorted.bam.bai)
    INSPECT_BAM_WORKLOAD(bam, inspector)
    inspected_bam = INSPECT_BAM_WORKLOAD.out.workload
    if (params.shard_mode == 'contig') {
        SPLIT_BAM_BY_CONTIG(inspected_bam, splitter)
        contig_bams = SPLIT_BAM_BY_CONTIG.out.shards.combine(SPLIT_BAM_BY_CONTIG.out.manifest).flatMap { meta, shard_dir, manifest ->
            manifest.readLines().drop(1).findAll { it.trim() }.collect { line ->
                def fields = line.split('\\t', -1)
                def contig = fields[0]
                def resource_class = fields[3]
                def mapped_reads = fields[2].toLong()
                def bam_file = file("${shard_dir}/${fields[4]}")
                tuple(meta, contig, resource_class, mapped_reads, bam_file, file("${bam_file}.bai"))
            }
        }
        TAMA_COLLAPSE(contig_bams, reference)
        VALIDATE_TAMA_OUTPUT(TAMA_COLLAPSE.out.bed, validator)
        accession_beds = VALIDATE_TAMA_OUTPUT.out.bed.map { meta, _shard, bed -> tuple(meta.id, bed) }
            .groupTuple()
            .map { accession, beds -> tuple(accession, beds.sort { left, right -> left.name <=> right.name }) }
        TAMA_MERGE_ACCESSION(accession_beds, merge_filelist_builder)
        final_beds = TAMA_MERGE_ACCESSION.out.bed.collect()
        version_ch = INSPECT_BAM_WORKLOAD.out.versions
            .mix(SPLIT_BAM_BY_CONTIG.out.versions)
            .mix(TAMA_COLLAPSE.out.versions)
            .mix(VALIDATE_TAMA_OUTPUT.out.versions)
            .mix(TAMA_MERGE_ACCESSION.out.versions)
    } else {
        whole_bams = inspected_bam.map { meta, bam_file, bai, workload ->
            def mapped_reads = workload.readLines().drop(1).findAll { it.trim() }.collect { it.split('\\t', -1)[2].toLong() }.sum()
            def resource_class = mapped_reads >= params.shard_contig_reads ? 'large' : 'small'
            tuple(meta, 'whole', resource_class, mapped_reads, bam_file, bai)
        }
        TAMA_COLLAPSE(whole_bams, reference)
        VALIDATE_TAMA_OUTPUT(TAMA_COLLAPSE.out.bed, validator)
        final_beds = VALIDATE_TAMA_OUTPUT.out.bed.map { _meta, _shard, bed -> bed }.collect()
        version_ch = INSPECT_BAM_WORKLOAD.out.versions
            .mix(TAMA_COLLAPSE.out.versions)
            .mix(VALIDATE_TAMA_OUTPUT.out.versions)
    }

    emit:
    beds = final_beds
    collapse_reports = INSPECT_BAM_WORKLOAD.out.workload
        .map { _meta, _bam, _bai, workload -> workload }
        .mix(TAMA_COLLAPSE.out.read)
        .mix(VALIDATE_TAMA_OUTPUT.out.report)
    versions = version_ch
}
