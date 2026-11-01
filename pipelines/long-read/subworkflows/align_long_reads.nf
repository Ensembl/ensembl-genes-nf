include { BUILD_MINIMAP_INDEX; ALIGN_FASTQ; SORT_FASTQ_ALIGNMENT; SORT_INPUT_BAM; ALIGNMENT_QC; SAMTOOLS_VIEW_ALIGNMENT; EXTRACT_ALIGNMENT_FEATURES; VALIDATE_ALIGNMENT_REFERENCE } from '../modules/alignment.nf'

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
    fastq_inputs = input_runs.filter { meta, _reads -> meta.read_representation == 'fastq' }
    bam_inputs = input_runs.filter { meta, _bam -> meta.read_representation in ['bam', 'ccs'] }
    ALIGN_FASTQ(fastq_inputs, reference, BUILD_MINIMAP_INDEX.out.index, minimap2_preset, secondary_mode)
    SORT_FASTQ_ALIGNMENT(ALIGN_FASTQ.out.sam, reference, minimap2_preset, secondary_mode)
    SORT_INPUT_BAM(bam_inputs, reference)
    raw_alignments = SORT_FASTQ_ALIGNMENT.out.alignment.mix(SORT_INPUT_BAM.out.alignment)
    VALIDATE_ALIGNMENT_REFERENCE(raw_alignments, BUILD_MINIMAP_INDEX.out.dictionary)
    alignments = VALIDATE_ALIGNMENT_REFERENCE.out.alignment
    ALIGNMENT_QC(alignments)
    SAMTOOLS_VIEW_ALIGNMENT(alignments)
    EXTRACT_ALIGNMENT_FEATURES(SAMTOOLS_VIEW_ALIGNMENT.out.sam)

    emit:
    alignment = alignments
    qc = ALIGNMENT_QC.out.qc
    features = EXTRACT_ALIGNMENT_FEATURES.out.features
    alignment_manifest = SORT_FASTQ_ALIGNMENT.out.manifest.mix(SORT_INPUT_BAM.out.manifest)
    versions = BUILD_MINIMAP_INDEX.out.versions.mix(ALIGN_FASTQ.out.versions).mix(SORT_FASTQ_ALIGNMENT.out.versions).mix(SORT_INPUT_BAM.out.versions).mix(VALIDATE_ALIGNMENT_REFERENCE.out.versions).mix(ALIGNMENT_QC.out.versions).mix(SAMTOOLS_VIEW_ALIGNMENT.out.versions).mix(EXTRACT_ALIGNMENT_FEATURES.out.versions)
    reference_dictionary = BUILD_MINIMAP_INDEX.out.dictionary
    reference_identity = BUILD_MINIMAP_INDEX.out.identity
    reference_validation = VALIDATE_ALIGNMENT_REFERENCE.out.validation
}
