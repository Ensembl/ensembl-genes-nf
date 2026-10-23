nextflow.enable.dsl = 2

include { BUILD_MINIMAP2_INDEX } from '../modules/build_minimap_index.nf'
include { MINIMAP2_ALIGN } from '../modules/minimap2_align.nf'
include { ALIGNMENT_QC } from '../modules/alignment_qc.nf'
include { ALIGNMENT_COUNTS } from '../modules/alignment_counts.nf'
include { ALIGNMENT_IDENTITY } from '../modules/alignment_identity.nf'
include { PER_BASE_COVERAGE } from '../modules/per_base_coverage.nf'
include { EXTRACT_INTRON_CHAINS } from '../modules/extract_intron_chains.nf'
include { SUMMARISE_SAMPLE_QC } from '../modules/summarise_sample_qc.nf'

workflow ALIGN_LONG_READS {
    take:
    reads
    pre_qc
    reference

    main:
    // reads: tuple val(meta), path(reads.fastq.gz)
    // bam: tuple val(meta), path(sorted.bam), path(sorted.bam.bai)
    build_versions = channel.empty()
    if (params.minimap_index) {
        minimap_index = channel.value(file(params.minimap_index, checkIfExists: true))
    } else {
        BUILD_MINIMAP2_INDEX(reference)
        // The reference is a single value, so the single generated index is
        // also reusable by every alignment task when paired with the read stream.
        minimap_index = BUILD_MINIMAP2_INDEX.out.index
        build_versions = BUILD_MINIMAP2_INDEX.out.versions
    }
    MINIMAP2_ALIGN(reads, minimap_index)
    qc_stats = channel.empty()
    alignment_counts = channel.empty()
    alignment_identity = channel.empty()
    coverage = channel.empty()
    intron_chains = channel.empty()
    sample_qc = channel.empty()
    qc_versions = channel.empty()
    if (params.alignment_qc == true || params.alignment_qc?.toString()?.toLowerCase() == 'true') {
        ALIGNMENT_QC(MINIMAP2_ALIGN.out.bam)
        identity_script = channel.value(file("${projectDir}/bin/alignment_identity.py", checkIfExists: true))
        chain_script = channel.value(file("${projectDir}/bin/extract_intron_chains.py", checkIfExists: true))
        ALIGNMENT_COUNTS(MINIMAP2_ALIGN.out.bam)
        ALIGNMENT_IDENTITY(MINIMAP2_ALIGN.out.bam, identity_script)
        PER_BASE_COVERAGE(MINIMAP2_ALIGN.out.bam)
        EXTRACT_INTRON_CHAINS(MINIMAP2_ALIGN.out.bam, chain_script)
        qc_stats = ALIGNMENT_QC.out.stats
        alignment_counts = ALIGNMENT_COUNTS.out.stats
        alignment_identity = ALIGNMENT_IDENTITY.out.stats
        coverage = PER_BASE_COVERAGE.out.summary
        intron_chains = EXTRACT_INTRON_CHAINS.out.summary
        qc_versions = ALIGNMENT_QC.out.versions
            .mix(ALIGNMENT_COUNTS.out.versions)
            .mix(ALIGNMENT_IDENTITY.out.versions)
            .mix(PER_BASE_COVERAGE.out.versions)
            .mix(EXTRACT_INTRON_CHAINS.out.versions)
        summary_script = channel.value(file("${projectDir}/bin/summarise_sample_qc.py", checkIfExists: true))
        sample_qc_input = pre_qc
            .join(alignment_counts, by: 0)
            .join(alignment_identity, by: 0)
            .join(coverage, by: 0)
            .join(intron_chains, by: 0)
            .map { meta, pre, counts, identity, cov, chains -> tuple(meta, pre, counts, identity, cov, chains) }
        SUMMARISE_SAMPLE_QC(sample_qc_input, summary_script)
        sample_qc = SUMMARISE_SAMPLE_QC.out.summary
        qc_versions = qc_versions.mix(SUMMARISE_SAMPLE_QC.out.versions)
    }

    emit:
    bam = MINIMAP2_ALIGN.out.bam
    reports = pre_qc.mix(MINIMAP2_ALIGN.out.config)
    stats = qc_stats
    pre_qc = pre_qc
    sample_qc = sample_qc
    versions = build_versions.mix(MINIMAP2_ALIGN.out.versions).mix(qc_versions)
}
