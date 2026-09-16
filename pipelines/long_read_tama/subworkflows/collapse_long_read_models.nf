nextflow.enable.dsl = 2

include { SPLIT_BAM_BY_CONTIG } from '../modules/split_contigs.nf'
include { TAMA_COLLAPSE } from '../modules/tama_collapse.nf'
include { TAMA_MERGE as TAMA_MERGE_ACCESSION } from '../modules/tama_merge.nf'
include { VALIDATE_TAMA_OUTPUT } from '../modules/validate_tama_output.nf'

workflow COLLAPSE_LONG_READ_MODELS {
    take:
    bam
    reference

    main:
    // bam: tuple val(meta), path(sorted.bam), path(sorted.bam.bai)
    if (params.shard_mode == 'contig') {
        SPLIT_BAM_BY_CONTIG(bam)
        contig_bams = SPLIT_BAM_BY_CONTIG.out.shards.flatMap { meta, shard_dir ->
            shard_dir.listFiles().findAll { bam_file -> bam_file.name.endsWith('.bam') }.collect {
                bam_file -> tuple(meta, bam_file.baseName, bam_file)
            }
        }
        TAMA_COLLAPSE(contig_bams, reference)
        VALIDATE_TAMA_OUTPUT(TAMA_COLLAPSE.out.bed)
        accession_beds = VALIDATE_TAMA_OUTPUT.out.bed.map { meta, _shard, bed -> tuple(meta.id, bed) }
            .groupTuple()
            .map { accession, beds -> tuple(accession, beds.sort { left, right -> left.name <=> right.name }) }
        TAMA_MERGE_ACCESSION(accession_beds)
        final_beds = TAMA_MERGE_ACCESSION.out.bed.collect()
        version_ch = SPLIT_BAM_BY_CONTIG.out.versions
            .mix(TAMA_COLLAPSE.out.versions)
            .mix(VALIDATE_TAMA_OUTPUT.out.versions)
            .mix(TAMA_MERGE_ACCESSION.out.versions)
    } else {
        whole_bams = bam.map { meta, bam_file, _bai -> tuple(meta, 'whole', bam_file) }
        TAMA_COLLAPSE(whole_bams, reference)
        VALIDATE_TAMA_OUTPUT(TAMA_COLLAPSE.out.bed)
        final_beds = VALIDATE_TAMA_OUTPUT.out.bed.map { _meta, _shard, bed -> bed }.collect()
        version_ch = TAMA_COLLAPSE.out.versions.mix(VALIDATE_TAMA_OUTPUT.out.versions)
    }

    emit:
    beds = final_beds
    collapse_reports = TAMA_COLLAPSE.out.read.mix(VALIDATE_TAMA_OUTPUT.out.report)
    versions = version_ch
}
