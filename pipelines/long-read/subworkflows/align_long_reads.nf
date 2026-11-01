include { BUILD_MINIMAP_INDEX; BUILD_REFERENCE_IDENTITY; ALIGN_FASTQ; SORT_FASTQ_ALIGNMENT; SORT_INPUT_BAM; WRITE_ALIGNMENT_MANIFEST as WRITE_FASTQ_ALIGNMENT_MANIFEST; WRITE_ALIGNMENT_MANIFEST as WRITE_BAM_ALIGNMENT_MANIFEST; SAMTOOLS_IDXSTATS; ALIGNMENT_QC; SAMTOOLS_VIEW_ALIGNMENT; EXTRACT_ALIGNMENT_FEATURES; VALIDATE_ALIGNMENT_REFERENCE } from '../modules/alignment.nf'

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
    WRITE_FASTQ_ALIGNMENT_MANIFEST(SORT_FASTQ_ALIGNMENT.out.manifest_input, reference, minimap2_preset, secondary_mode)
    WRITE_BAM_ALIGNMENT_MANIFEST(SORT_INPUT_BAM.out.manifest_input, reference, '', false)
    raw_alignments = SORT_FASTQ_ALIGNMENT.out.alignment.mix(SORT_INPUT_BAM.out.alignment)
    SAMTOOLS_IDXSTATS(raw_alignments)
    VALIDATE_ALIGNMENT_REFERENCE(SAMTOOLS_IDXSTATS.out.idxstats, BUILD_REFERENCE_IDENTITY.out.dictionary)
    alignments = VALIDATE_ALIGNMENT_REFERENCE.out.alignment
    ALIGNMENT_QC(alignments)
    SAMTOOLS_VIEW_ALIGNMENT(alignments)
    EXTRACT_ALIGNMENT_FEATURES(SAMTOOLS_VIEW_ALIGNMENT.out.sam)

    emit:
    alignment = alignments
    qc = ALIGNMENT_QC.out.qc
    features = EXTRACT_ALIGNMENT_FEATURES.out.features
    alignment_manifest = WRITE_FASTQ_ALIGNMENT_MANIFEST.out.manifest.mix(WRITE_BAM_ALIGNMENT_MANIFEST.out.manifest)
    versions = BUILD_MINIMAP_INDEX.out.versions.mix(BUILD_REFERENCE_IDENTITY.out.versions).mix(ALIGN_FASTQ.out.versions).mix(SORT_FASTQ_ALIGNMENT.out.versions).mix(SORT_INPUT_BAM.out.versions).mix(WRITE_FASTQ_ALIGNMENT_MANIFEST.out.versions).mix(WRITE_BAM_ALIGNMENT_MANIFEST.out.versions).mix(SAMTOOLS_IDXSTATS.out.versions).mix(VALIDATE_ALIGNMENT_REFERENCE.out.versions).mix(ALIGNMENT_QC.out.versions).mix(SAMTOOLS_VIEW_ALIGNMENT.out.versions).mix(EXTRACT_ALIGNMENT_FEATURES.out.versions)
    reference_dictionary = BUILD_REFERENCE_IDENTITY.out.dictionary
    reference_identity = BUILD_REFERENCE_IDENTITY.out.identity
    reference_validation = VALIDATE_ALIGNMENT_REFERENCE.out.validation
}
