include { VALIDATE_INPUT_MANIFEST; FINALIZE_PREPARED_MANIFEST; VALIDATE_FASTQ; AUDIT_FASTQ; VALIDATE_BAM_INPUT; STAGE_VALIDATED_INPUT; STAGE_VALIDATED_FASTQ; DISCOVER_LONG_READ_DATA; RESOLVE_LONG_READ_METADATA; ACQUIRE_FASTQ; ACQUIRE_BAM; RUN_PBCCS; INDEX_PBCCS; BAM_TO_FASTQ } from '../modules/preparation.nf'

workflow PREPARE_LONG_READS {
    take:
    manifest
    input_runs
    mode
    cohort_id
    manifest_root
    ccs_mode
    preparation_stages
    discovery_query
    cache_dir

    main:
    allow_unresolved_inputs = preparation_stages.contains('acquire')
    VALIDATE_INPUT_MANIFEST(manifest, mode, cohort_id, manifest_root, ccs_mode, allow_unresolved_inputs, cache_dir)
    if (preparation_stages.contains('discover')) {
        DISCOVER_LONG_READ_DATA(discovery_query ?: manifest.toString())
    }
    if (preparation_stages.contains('resolve_metadata')) {
        RESOLVE_LONG_READ_METADATA(manifest)
    }

    acquired_fastq = channel.empty()
    acquired_bam = channel.empty()
    if (preparation_stages.contains('acquire')) {
        acquire_fastq_requests = input_runs.filter { item -> item[0].read_representation == 'fastq' }
            .map { meta, _source -> tuple(meta, meta.source_uri) }
        acquire_bam_requests = input_runs.filter { item -> item[0].read_representation in ['bam', 'ccs', 'raw_subreads'] }
            .map { meta, _source -> tuple(meta, meta.source_uri) }
        ACQUIRE_FASTQ(acquire_fastq_requests)
        ACQUIRE_BAM(acquire_bam_requests)
        acquired_fastq = ACQUIRE_FASTQ.out.reads
        acquired_bam = ACQUIRE_BAM.out.bam
    }

    fastq_inputs = preparation_stages.contains('acquire')
        ? acquired_fastq
        : input_runs.filter { item -> item[0].read_representation == 'fastq' }
    bam_inputs = preparation_stages.contains('acquire')
        ? acquired_bam
        : input_runs.filter { item -> item[0].read_representation in ['bam', 'ccs', 'raw_subreads'] }
    VALIDATE_BAM_INPUT(bam_inputs)

    bam_to_fastq_versions = channel.empty()
    if (preparation_stages.contains('bam_to_fastq')) {
        if (ccs_mode == 'pbccs') {
            error 'preparation_stages bam_to_fastq cannot be combined with ccs_mode pbccs'
        }
        BAM_TO_FASTQ(VALIDATE_BAM_INPUT.out.bam.map { meta, bam, _bai, _status -> tuple(meta, bam) })
        fastq_inputs = BAM_TO_FASTQ.out.fastq
        bam_to_fastq_versions = BAM_TO_FASTQ.out.versions
    }
    VALIDATE_FASTQ(fastq_inputs)
    AUDIT_FASTQ(VALIDATE_FASTQ.out.reads)
    STAGE_VALIDATED_FASTQ(AUDIT_FASTQ.out.reads)

    ccs_runs = VALIDATE_BAM_INPUT.out.bam
        .filter { meta, _bam, _bai, _status -> meta.read_representation == 'raw_subreads' && ccs_mode == 'pbccs' }
        .map { meta, bam, bai, _status -> tuple(meta, bam, bai) }
    RUN_PBCCS(ccs_runs)
    INDEX_PBCCS(RUN_PBCCS.out.ccs)
    validated_bam = preparation_stages.contains('bam_to_fastq')
        ? channel.empty()
        : VALIDATE_BAM_INPUT.out.bam.filter { meta, _bam, _bai, _status -> meta.read_representation == 'bam' }
            .map { meta, bam, _bai, _status -> tuple(meta, bam) }
    validated_ccs = preparation_stages.contains('bam_to_fastq')
        ? channel.empty()
        : ccs_mode == 'pbccs'
            ? INDEX_PBCCS.out.ccs.map { meta, bam, _pbi -> tuple(meta, bam) }
            : VALIDATE_BAM_INPUT.out.bam.filter { meta, _bam, _bai, _status -> meta.read_representation == 'ccs' }
                .map { meta, bam, _bai, _status -> tuple(meta, bam) }
    STAGE_VALIDATED_INPUT(validated_bam.mix(validated_ccs))

    prepared_metadata = channel.empty()
    if (preparation_stages.contains('acquire') || preparation_stages.contains('bam_to_fastq')) {
        prepared_metadata = preparation_stages.contains('bam_to_fastq')
            ? STAGE_VALIDATED_FASTQ.out.prepared
            : STAGE_VALIDATED_FASTQ.out.prepared.mix(VALIDATE_BAM_INPUT.out.prepared)
        prepared_metadata = prepared_metadata.collectFile(name: 'prepared_inputs.tsv', newLine: true)
        FINALIZE_PREPARED_MANIFEST(VALIDATE_INPUT_MANIFEST.out.normalized, prepared_metadata)
    }

    versions = VALIDATE_BAM_INPUT.out.versions
        .mix(RUN_PBCCS.out.versions)
        .mix(bam_to_fastq_versions)

    emit:
    reads = STAGE_VALIDATED_FASTQ.out.reads.map { meta, reads -> tuple(meta, reads) }
        .mix(STAGE_VALIDATED_INPUT.out.reads)
    normalized = (preparation_stages.contains('acquire') || preparation_stages.contains('bam_to_fastq'))
        ? FINALIZE_PREPARED_MANIFEST.out.normalized
        : VALIDATE_INPUT_MANIFEST.out.normalized
    versions
}
