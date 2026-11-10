include { BUILD_MINIMAP_INDEX; BUILD_REFERENCE_IDENTITY; ALIGN_FASTQ; SORT_FASTQ_ALIGNMENT; SORT_INPUT_BAM; INDEX_SORTED_ALIGNMENT as INDEX_FASTQ_ALIGNMENT; INDEX_SORTED_ALIGNMENT as INDEX_BAM_ALIGNMENT; VALIDATE_SORTED_ALIGNMENT as VALIDATE_FASTQ_ALIGNMENT; VALIDATE_SORTED_ALIGNMENT as VALIDATE_BAM_ALIGNMENT; WRITE_ALIGNMENT_MANIFEST as WRITE_FASTQ_ALIGNMENT_MANIFEST; WRITE_ALIGNMENT_MANIFEST as WRITE_BAM_ALIGNMENT_MANIFEST; SAMTOOLS_IDXSTATS; SAMTOOLS_FLAGSTAT; SAMTOOLS_QC_SUMMARY; SAMTOOLS_VIEW_ALIGNMENT; EXTRACT_ALIGNMENT_FEATURES; VALIDATE_ALIGNMENT_REFERENCE } from '../modules/alignment.nf'

workflow ALIGN_LONG_READS {
    take:
    input_runs
    reference
    _input_mode
    minimap2_preset
    secondary_mode
    _run_qc

    main:
    BUILD_MINIMAP_INDEX(reference)
    BUILD_REFERENCE_IDENTITY(reference)
    fastq_inputs = input_runs.filter { meta, _reads -> meta.read_representation == 'fastq' }
    bam_inputs = input_runs.filter { meta, _bam -> meta.read_representation in ['bam', 'ccs'] }
    ALIGN_FASTQ(fastq_inputs, reference, BUILD_MINIMAP_INDEX.out.index, minimap2_preset, secondary_mode)
    SORT_FASTQ_ALIGNMENT(ALIGN_FASTQ.out.sam, reference, minimap2_preset, secondary_mode)
    SORT_INPUT_BAM(bam_inputs, reference)
    INDEX_FASTQ_ALIGNMENT(SORT_FASTQ_ALIGNMENT.out.raw)
    INDEX_BAM_ALIGNMENT(SORT_INPUT_BAM.out.raw)
    VALIDATE_FASTQ_ALIGNMENT(INDEX_FASTQ_ALIGNMENT.out.indexed)
    VALIDATE_BAM_ALIGNMENT(INDEX_BAM_ALIGNMENT.out.indexed)
    WRITE_FASTQ_ALIGNMENT_MANIFEST(VALIDATE_FASTQ_ALIGNMENT.out.manifest_input, reference, minimap2_preset, secondary_mode)
    WRITE_BAM_ALIGNMENT_MANIFEST(VALIDATE_BAM_ALIGNMENT.out.manifest_input, reference, '', false)
    raw_alignments = VALIDATE_FASTQ_ALIGNMENT.out.alignment.mix(VALIDATE_BAM_ALIGNMENT.out.alignment)
    SAMTOOLS_IDXSTATS(raw_alignments)
    VALIDATE_ALIGNMENT_REFERENCE(SAMTOOLS_IDXSTATS.out.idxstats, BUILD_REFERENCE_IDENTITY.out.dictionary)
    alignments = VALIDATE_ALIGNMENT_REFERENCE.out.alignment
    SAMTOOLS_FLAGSTAT(alignments)
    qc_inputs = alignments
        .join(SAMTOOLS_IDXSTATS.out.idxstats.map { meta, _bam, _bai, idxstats -> tuple(meta, idxstats) }, by: 0)
        .map { meta, bam, bai, idxstats -> tuple(meta, bam, bai, idxstats) }
    SAMTOOLS_QC_SUMMARY(qc_inputs)
    SAMTOOLS_VIEW_ALIGNMENT(alignments)
    EXTRACT_ALIGNMENT_FEATURES(SAMTOOLS_VIEW_ALIGNMENT.out.sam)

    emit:
    alignment = alignments
    qc = SAMTOOLS_FLAGSTAT.out.flagstat
        .join(SAMTOOLS_IDXSTATS.out.idxstats.map { meta, _bam, _bai, idxstats -> tuple(meta, idxstats) }, by: 0)
        .join(SAMTOOLS_QC_SUMMARY.out.qc, by: 0)
        .map { meta, flagstat, idxstats, quickcheck, alignment_summary -> tuple(meta, flagstat, idxstats, quickcheck, alignment_summary) }
    features = EXTRACT_ALIGNMENT_FEATURES.out.features
    direct_support = EXTRACT_ALIGNMENT_FEATURES.out.direct_support
    alignment_manifest = WRITE_FASTQ_ALIGNMENT_MANIFEST.out.manifest.mix(WRITE_BAM_ALIGNMENT_MANIFEST.out.manifest)
    versions = BUILD_MINIMAP_INDEX.out.versions.mix(ALIGN_FASTQ.out.versions).mix(SORT_FASTQ_ALIGNMENT.out.versions).mix(SORT_INPUT_BAM.out.versions).mix(INDEX_FASTQ_ALIGNMENT.out.versions).mix(INDEX_BAM_ALIGNMENT.out.versions).mix(VALIDATE_FASTQ_ALIGNMENT.out.versions).mix(VALIDATE_BAM_ALIGNMENT.out.versions).mix(SAMTOOLS_IDXSTATS.out.versions).mix(SAMTOOLS_FLAGSTAT.out.versions).mix(SAMTOOLS_QC_SUMMARY.out.versions).mix(SAMTOOLS_VIEW_ALIGNMENT.out.versions)
    reference_dictionary = BUILD_REFERENCE_IDENTITY.out.dictionary
    reference_identity = BUILD_REFERENCE_IDENTITY.out.identity
    reference_validation = VALIDATE_ALIGNMENT_REFERENCE.out.validation
}
