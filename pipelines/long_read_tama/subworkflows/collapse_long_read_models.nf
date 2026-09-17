nextflow.enable.dsl = 2

include { SPLIT_BAM_BY_CONTIG } from '../modules/split_contigs.nf'
include { INSPECT_BAM_WORKLOAD } from '../modules/inspect_bam_workload.nf'
include { TAMA_COLLAPSE } from '../modules/tama_collapse.nf'
include { STRINGTIE2_COLLAPSE } from '../modules/stringtie2_collapse.nf'
include { STRINGTIE3_COLLAPSE } from '../modules/stringtie3_collapse.nf'
include { TMERGE_COLLAPSE } from '../modules/tmerge_collapse.nf'
include { TAMA_MERGE as TAMA_MERGE_ACCESSION } from '../modules/tama_merge.nf'
include { TMERGE as TMERGE_ACCESSION } from '../modules/tmerge.nf'
include { VALIDATE_TAMA_OUTPUT } from '../modules/validate_tama_output.nf'

workflow COLLAPSE_LONG_READ_MODELS {
    take:
    bam
    reference

    main:
    // bam: tuple val(meta), path(sorted.bam), path(sorted.bam.bai)
    INSPECT_BAM_WORKLOAD(bam)
    inspected_bam = INSPECT_BAM_WORKLOAD.out.workload
    if (params.shard_mode == 'contig') {
        SPLIT_BAM_BY_CONTIG(inspected_bam)
        contig_bams = SPLIT_BAM_BY_CONTIG.out.shards.flatMap { meta, shard_dir, manifest ->
            manifest.readLines().drop(1).findAll { it.trim() }.collect { line ->
                def fields = line.split('\\t', -1)
                def contig = fields[0]
                def resource_class = fields[3]
                def mapped_reads = fields[2].toLong()
                def bam_file = file("${shard_dir}/${fields[4]}")
                tuple(meta, contig, resource_class, mapped_reads, bam_file, file("${bam_file}.bai"))
            }
        }
        run_tama = params.model_backend in ['tama', 'all']
        run_stringtie2 = params.model_backend in ['stringtie2', 'all']
        run_stringtie3 = params.model_backend in ['stringtie3', 'all']
        run_tmerge = params.model_backend in ['tmerge', 'all']

        if (run_tama) {
            TAMA_COLLAPSE(contig_bams, reference)
            VALIDATE_TAMA_OUTPUT(TAMA_COLLAPSE.out.bed)
        }
        if (run_stringtie2) {
            STRINGTIE2_COLLAPSE(contig_bams)
        }
        if (run_stringtie3) {
            STRINGTIE3_COLLAPSE(contig_bams)
        }
        if (run_tmerge) {
            TMERGE_COLLAPSE(contig_bams)
        }

        // In comparison mode TAMA remains the canonical downstream model set;
        // the other backend outputs are published independently for review.
        selected_beds = run_tama ? VALIDATE_TAMA_OUTPUT.out.bed :
            (run_stringtie2 ? STRINGTIE2_COLLAPSE.out.bed :
            (run_stringtie3 ? STRINGTIE3_COLLAPSE.out.bed : TMERGE_COLLAPSE.out.bed))
        accession_beds = selected_beds.map { meta, _shard, bed -> tuple(meta.id, bed) }
            .groupTuple()
            .map { accession, beds -> tuple(accession, beds.sort { left, right -> left.name <=> right.name }) }
        if (params.merge_tool == 'tmerge' && run_tama) {
            TMERGE_ACCESSION(accession_beds)
            final_beds = TMERGE_ACCESSION.out.bed.collect()
        } else {
            TAMA_MERGE_ACCESSION(accession_beds)
            final_beds = TAMA_MERGE_ACCESSION.out.bed.collect()
        }
        version_ch = INSPECT_BAM_WORKLOAD.out.versions
            .mix(SPLIT_BAM_BY_CONTIG.out.versions)
        if (run_tama) {
            version_ch = version_ch.mix(TAMA_COLLAPSE.out.versions).mix(VALIDATE_TAMA_OUTPUT.out.versions)
        }
        if (run_stringtie2) { version_ch = version_ch.mix(STRINGTIE2_COLLAPSE.out.versions) }
        if (run_stringtie3) { version_ch = version_ch.mix(STRINGTIE3_COLLAPSE.out.versions) }
        if (run_tmerge) { version_ch = version_ch.mix(TMERGE_COLLAPSE.out.versions) }
        version_ch = version_ch.mix(params.merge_tool == 'tmerge' && run_tama ? TMERGE_ACCESSION.out.versions : TAMA_MERGE_ACCESSION.out.versions)
    } else {
        whole_bams = inspected_bam.map { meta, bam_file, bai, workload ->
            def mapped_reads = workload.readLines().drop(1).findAll { it.trim() }.collect { it.split('\\t', -1)[2].toLong() }.sum()
            def resource_class = mapped_reads >= params.shard_contig_reads ? 'large' : 'small'
            tuple(meta, 'whole', resource_class, mapped_reads, bam_file, bai)
        }
        run_tama = params.model_backend in ['tama', 'all']
        run_stringtie2 = params.model_backend in ['stringtie2', 'all']
        run_stringtie3 = params.model_backend in ['stringtie3', 'all']
        run_tmerge = params.model_backend in ['tmerge', 'all']
        if (run_tama) { TAMA_COLLAPSE(whole_bams, reference); VALIDATE_TAMA_OUTPUT(TAMA_COLLAPSE.out.bed) }
        if (run_stringtie2) { STRINGTIE2_COLLAPSE(whole_bams) }
        if (run_stringtie3) { STRINGTIE3_COLLAPSE(whole_bams) }
        if (run_tmerge) { TMERGE_COLLAPSE(whole_bams) }
        selected_beds = run_tama ? VALIDATE_TAMA_OUTPUT.out.bed :
            (run_stringtie2 ? STRINGTIE2_COLLAPSE.out.bed :
            (run_stringtie3 ? STRINGTIE3_COLLAPSE.out.bed : TMERGE_COLLAPSE.out.bed))
        final_beds = selected_beds.map { _meta, _shard, bed -> bed }.collect()
        version_ch = INSPECT_BAM_WORKLOAD.out.versions
        if (run_tama) { version_ch = version_ch.mix(TAMA_COLLAPSE.out.versions).mix(VALIDATE_TAMA_OUTPUT.out.versions) }
        if (run_stringtie2) { version_ch = version_ch.mix(STRINGTIE2_COLLAPSE.out.versions) }
        if (run_stringtie3) { version_ch = version_ch.mix(STRINGTIE3_COLLAPSE.out.versions) }
        if (run_tmerge) { version_ch = version_ch.mix(TMERGE_COLLAPSE.out.versions) }
    }

    collapse_reports = INSPECT_BAM_WORKLOAD.out.workload.map { _meta, _bam, _bai, workload -> workload }
    if (params.model_backend in ['tama', 'all']) {
        collapse_reports = collapse_reports.mix(TAMA_COLLAPSE.out.read).mix(TAMA_COLLAPSE.out.status).mix(TAMA_COLLAPSE.out.stderr).mix(VALIDATE_TAMA_OUTPUT.out.report)
    }

    emit:
    beds = final_beds
    collapse_reports
    versions = version_ch
}
