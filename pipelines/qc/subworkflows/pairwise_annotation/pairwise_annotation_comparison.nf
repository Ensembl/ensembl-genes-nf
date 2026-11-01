include { VALIDATE_PAIR_MANIFEST } from '../../modules/pairwise_annotation/validate_pair_manifest.nf'
include { STAGE_ANNOTATION_PAIR } from '../../modules/pairwise_annotation/stage_annotation_pair.nf'
include { NORMALISE_TRANSCRIPT_MODELS } from '../../modules/pairwise_annotation/normalise_transcript_models.nf'
include { ASSIGN_TRANSCRIPT_LOCI } from '../../modules/pairwise_annotation/assign_transcript_loci.nf'
include { AUDIT_LOCUS_ASSIGNMENT } from '../../modules/pairwise_annotation/audit_locus_assignment.nf'
include { RENDER_SYNTHETIC_ANNOTATION } from '../../modules/pairwise_annotation/render_synthetic_annotation.nf'
include { AUDIT_ANNOTATION_CONVERSION } from '../../modules/pairwise_annotation/audit_annotation_conversion.nf'
include { SUMMARISE_MODEL_COUNTS } from '../../modules/pairwise_annotation/summarise_model_counts.nf'
include { COMPARE_INTRON_CHAINS } from '../../modules/pairwise_annotation/compare_intron_chains.nf'
include { COMPARE_TRANSCRIPT_ENDS } from '../../modules/pairwise_annotation/compare_transcript_ends.nf'
include { COMPARE_UTR_STRUCTURE } from '../../modules/pairwise_annotation/compare_utr_structure.nf'
include { SUMMARISE_STRUCTURAL_COMPARISON } from '../../modules/pairwise_annotation/summarise_structural_comparison.nf'
include { METRIC_STATUS } from '../../modules/pairwise_annotation/metric_status.nf'
include { ASSEMBLE_PAIR_REPORT } from '../../modules/pairwise_annotation/assemble_pair_report.nf'
include { COLLECT_SOFTWARE_VERSIONS } from '../../modules/pairwise_annotation/collect_software_versions.nf'
include { ANNOTATE_MODEL_EVIDENCE } from '../../modules/pairwise_annotation/annotate_model_evidence.nf'
include { PREPARE_GFFCOMPARE_INPUTS } from '../../modules/pairwise_annotation/prepare_gffcompare_inputs.nf'
include { GFFCOMPARE } from '../../modules/pairwise_annotation/gffcompare.nf'
include { PARSE_GFFCOMPARE_TMAP } from '../../modules/pairwise_annotation/parse_gffcompare_tmap.nf'
include { SUMMARISE_GFFCOMPARE_MODELS } from '../../modules/pairwise_annotation/summarise_gffcompare_models.nf'
include { PAIR_GENES } from '../../modules/pairwise_annotation/pair_genes.nf'
include { PAIRWISE_TRANSCRIPT_CONCORDANCE } from '../../modules/pairwise_annotation/transcript_concordance.nf'
include { PAIRWISE_CODING_INTEGRITY } from '../../modules/pairwise_annotation/coding_integrity.nf'
include { PAIRWISE_GENE_PRESENCE } from '../../modules/pairwise_annotation/gene_presence.nf'
include { PAIRWISE_MULTI_MAPPING } from '../../modules/pairwise_annotation/multi_mapping.nf'
include { PAIRWISE_GFF_FEATURE_METRICS } from '../../modules/pairwise_annotation/gff_feature_metrics.nf'
include { PAIRWISE_COUNT_GFF_TRANSCRIPTS } from '../../modules/pairwise_annotation/count_gff_transcripts.nf'

workflow PAIRWISE_ANNOTATION_COMPARISON {
    take:
    pairwise_samples_ch

    main:
    // Public input: tuple(meta, source_a, source_b, assembly_report,
    // evidence_a, evidence_b, reference_annotation).
    VALIDATE_PAIR_MANIFEST(pairwise_samples_ch)
    validated = VALIDATE_PAIR_MANIFEST.out.validated
    STAGE_ANNOTATION_PAIR(validated.map { meta, a, b, report, evidence_a, evidence_b, reference, _validated -> tuple(meta, a, b, report, evidence_a, evidence_b, reference) }, params.pairwise_cache_dir)
    staged = STAGE_ANNOTATION_PAIR.out.staged

    models_input = staged.map { meta, a, b, _report, _evidence_a, _evidence_b, _reference -> tuple(meta, a, b) }
    NORMALISE_TRANSCRIPT_MODELS(models_input)
    raw_models = NORMALISE_TRANSCRIPT_MODELS.out.models

    // The reference is staged as a value for every pair; an absent reference
    // is an empty file, so the assignment process always has an explicit mode.
    raw_keyed = raw_models.map { meta, a_models, b_models -> tuple(meta.pair_id, meta, a_models, b_models) }
    reference_keyed = staged.map { meta, _a, _b, _report, _ea, _eb, reference -> tuple(meta.pair_id, meta, reference) }
    models_with_reference = raw_keyed.join(reference_keyed).map { _pair_id, models_meta, a_models, b_models, _staged_meta, reference -> tuple(models_meta, a_models, b_models, reference) }
    ASSIGN_TRANSCRIPT_LOCI(models_with_reference)
    assigned = ASSIGN_TRANSCRIPT_LOCI.out.loci

    canonical_models = assigned.map { meta, a_models, b_models, _a_loci, _b_loci, _assignment_summary -> tuple(meta, a_models, b_models) }
    AUDIT_LOCUS_ASSIGNMENT(assigned.map { meta, a_models, b_models, _a_loci, _b_loci, assignment_summary -> tuple(meta, a_models, b_models, assignment_summary) })
    SUMMARISE_MODEL_COUNTS(canonical_models)

    RENDER_SYNTHETIC_ANNOTATION(assigned.flatMap { meta, a_models, b_models, _a_loci, _b_loci, _summary -> [tuple(meta, 'source_a', a_models), tuple(meta, 'source_b', b_models)] })
    AUDIT_ANNOTATION_CONVERSION(assigned.map { meta, a_models, b_models, _a_loci, _b_loci, _summary -> tuple(meta, a_models, b_models) })

    // Synthetic GFFs are valid inputs for structural feature/count metrics.
    // They must not be routed into biological gene, RBH, or CDS metrics.
    synthetic_a = RENDER_SYNTHETIC_ANNOTATION.out.annotation
        .filter { _meta, source_label, _gff -> source_label == 'source_a' }
        .map { meta, _source_label, gff -> tuple(meta.pair_id, meta, gff) }
    synthetic_b = RENDER_SYNTHETIC_ANNOTATION.out.annotation
        .filter { _meta, source_label, _gff -> source_label == 'source_b' }
        .map { meta, _source_label, gff -> tuple(meta.pair_id, meta, gff) }
    synthetic_pairs = synthetic_a.join(synthetic_b)
        .map { _pair_id, meta_a, source_a_gff, _meta_b, source_b_gff -> tuple(meta_a, source_a_gff, source_b_gff) }

    gff_feature_versions = channel.empty()
    transcript_count_versions = channel.empty()
    if (params.run_pairwise_gff_feature_metrics) {
        PAIRWISE_GFF_FEATURE_METRICS(synthetic_pairs)
        gff_feature_versions = PAIRWISE_GFF_FEATURE_METRICS.out.versions
    }
    if (params.run_pairwise_transcript_counts) {
        PAIRWISE_COUNT_GFF_TRANSCRIPTS(synthetic_pairs.flatMap { meta, source_a_gff, source_b_gff -> [
            tuple(meta, 'source_a', source_a_gff), tuple(meta, 'source_b', source_b_gff)
        ] })
        transcript_count_versions = PAIRWISE_COUNT_GFF_TRANSCRIPTS.out.versions
    }

    // Evidence joins to canonical transcript IDs and never to run-level files.
    canonical_keyed = canonical_models.map { meta, a_models, b_models -> tuple(meta.pair_id, meta, a_models, b_models) }
    evidence_keyed = staged.map { meta, _a, _b, _report, evidence_a, evidence_b, _reference -> tuple(meta.pair_id, meta, evidence_a, evidence_b) }
    canonical_keyed.join(evidence_keyed)
        .filter { _pair_id, joined_meta, _a_models, _b_models, _staged_meta, _evidence_a, _evidence_b -> joined_meta.evidence_a && joined_meta.evidence_b }
        .flatMap { _pair_id, meta, a_models, b_models, _staged_meta, evidence_a, evidence_b -> [tuple(meta, 'source_a', a_models, evidence_a), tuple(meta, 'source_b', b_models, evidence_b)] }
        .set { evidence_inputs }
    ANNOTATE_MODEL_EVIDENCE(evidence_inputs)

    COMPARE_INTRON_CHAINS(canonical_models)
    COMPARE_TRANSCRIPT_ENDS(canonical_models)
    COMPARE_UTR_STRUCTURE(canonical_models)

    intron_keyed = COMPARE_INTRON_CHAINS.out.results.map { meta, file -> tuple(meta.pair_id, meta, file) }
    ends_keyed = COMPARE_TRANSCRIPT_ENDS.out.results.map { meta, file -> tuple(meta.pair_id, meta, file) }
    utr_keyed = COMPARE_UTR_STRUCTURE.out.results.map { meta, file -> tuple(meta.pair_id, meta, file) }
    intron_keyed.join(ends_keyed).join(utr_keyed).map { _pair_id, meta, intron, _meta_ends, ends, _meta_utr, utr -> tuple(meta, intron, ends, utr) }.set { structural_summary_input }
    SUMMARISE_STRUCTURAL_COMPARISON(structural_summary_input)
    ASSEMBLE_PAIR_REPORT(SUMMARISE_STRUCTURAL_COMPARISON.out.summary.map { meta, summary -> tuple(meta, summary) })

    gffcompare_summary = channel.empty()
    gffcompare_versions = channel.empty()
    if (params.run_pairwise_gffcompare) {
        PREPARE_GFFCOMPARE_INPUTS(canonical_models)
        gffcompare_inputs = PREPARE_GFFCOMPARE_INPUTS.out.gtfs.flatMap { meta, source_a_gtf, source_b_gtf -> [
            tuple(meta, 'source_a_to_source_b', source_a_gtf, source_b_gtf),
            tuple(meta, 'source_b_to_source_a', source_b_gtf, source_a_gtf)
        ] }
        GFFCOMPARE(gffcompare_inputs)
        PARSE_GFFCOMPARE_TMAP(GFFCOMPARE.out.tmap)
        gffcompare_classes = PARSE_GFFCOMPARE_TMAP.out.classes.map { meta, _direction, classes -> tuple(meta.pair_id, meta, classes) }.groupTuple().map { _pair_id, metas, classes -> tuple(metas[0], classes) }
        SUMMARISE_GFFCOMPARE_MODELS(gffcompare_classes)
        gffcompare_summary = SUMMARISE_GFFCOMPARE_MODELS.out.summary
        gffcompare_versions = PREPARE_GFFCOMPARE_INPUTS.out.versions.mix(GFFCOMPARE.out.versions).mix(PARSE_GFFCOMPARE_TMAP.out.versions).mix(SUMMARISE_GFFCOMPARE_MODELS.out.versions)
    }

    status_input = pairwise_samples_ch.map { meta, _a, _b, _report, _ea, _eb, _reference ->
        def has_bed12 = meta.format_a == 'bed12' || meta.format_b == 'bed12'
        def structural = [
            ['transcript_counts', 'enabled', 'canonical transcript models'],
            ['exon_counts', 'enabled', 'canonical transcript models'],
            ['intron_chain', 'enabled', 'union-derived canonical loci'],
            ['transcript_ends', 'enabled', 'canonical transcript models'],
            ['utr_structure', 'enabled', 'structural UTR fields where available'],
            ['locus_counts', 'enabled', 'source and pair locus summaries'],
            ['gffcompare', params.run_pairwise_gffcompare ? 'enabled' : 'not_applicable', params.run_pairwise_gffcompare ? 'requested' : 'disabled by parameter']
        ]
        def biological = [
            ['gene_presence', has_bed12 ? 'not_applicable' : 'enabled', has_bed12 ? 'BED12 has no biological gene contract' : 'biological annotation input'],
            ['gene_rbh_pairing', has_bed12 ? 'not_applicable' : 'enabled', has_bed12 ? 'BED12 has no biological gene contract' : 'biological annotation input'],
            ['cds_coding_integrity', has_bed12 ? 'not_applicable' : 'missing_input', has_bed12 ? 'BED12 has no CDS contract' : 'CDS evidence must be supplied explicitly'],
            ['protein_validation', 'not_applicable', 'protein evidence is outside this component']
        ]
        tuple(meta, structural + biological)
    }
    METRIC_STATUS(status_input)

    pairwise_versions = VALIDATE_PAIR_MANIFEST.out.versions
        .mix(STAGE_ANNOTATION_PAIR.out.versions)
        .mix(NORMALISE_TRANSCRIPT_MODELS.out.versions)
        .mix(ASSIGN_TRANSCRIPT_LOCI.out.versions)
        .mix(AUDIT_LOCUS_ASSIGNMENT.out.versions)
        .mix(RENDER_SYNTHETIC_ANNOTATION.out.versions)
        .mix(AUDIT_ANNOTATION_CONVERSION.out.versions)
        .mix(SUMMARISE_MODEL_COUNTS.out.versions)
        .mix(COMPARE_INTRON_CHAINS.out.versions)
        .mix(COMPARE_TRANSCRIPT_ENDS.out.versions)
        .mix(COMPARE_UTR_STRUCTURE.out.versions)
        .mix(SUMMARISE_STRUCTURAL_COMPARISON.out.versions)
        .mix(METRIC_STATUS.out.versions)
        .mix(ASSEMBLE_PAIR_REPORT.out.versions)
        .mix(ANNOTATE_MODEL_EVIDENCE.out.versions)
        .mix(gff_feature_versions)
        .mix(transcript_count_versions)
        .mix(gffcompare_versions)
        .collect()
    COLLECT_SOFTWARE_VERSIONS(pairwise_versions)

    // Legacy biological metrics receive the original annotations only when
    // both inputs are annotation formats with biological feature contracts.
    legacy = staged.filter { meta, _a, _b, _report, _ea, _eb, _reference -> meta.format_a != 'bed12' && meta.format_b != 'bed12' }
    if (params.run_pairwise_legacy_metrics) {
        PAIR_GENES(legacy.map { meta, a, b, report, _ea, _eb, _reference -> tuple(meta, a, b, report) })
        rbh_inputs = legacy.map { meta, a, b, _report, _ea, _eb, _reference -> tuple(meta, a, b) }.join(PAIR_GENES.out.rbh)
        if (params.run_pairwise_transcript_concordance) PAIRWISE_TRANSCRIPT_CONCORDANCE(rbh_inputs)
        if (params.run_pairwise_coding_integrity) PAIRWISE_CODING_INTEGRITY(rbh_inputs)
        if (params.run_pairwise_gene_presence) PAIRWISE_GENE_PRESENCE(legacy.map { meta, a, b, _report, _ea, _eb, _reference -> tuple(meta, a, b) }, '')
        if (params.run_pairwise_multi_mapping) PAIRWISE_MULTI_MAPPING(PAIR_GENES.out.all_pairs)
    }

    emit:
    staged_inputs = staged
    raw_models = raw_models
    canonical_models = canonical_models
    locus_assignments = assigned
    locus_audit = AUDIT_LOCUS_ASSIGNMENT.out.audit
    conversion_audits = AUDIT_ANNOTATION_CONVERSION.out.audit
    id_maps = AUDIT_ANNOTATION_CONVERSION.out.id_map
    synthetic_annotations = RENDER_SYNTHETIC_ANNOTATION.out.annotation
    model_counts = SUMMARISE_MODEL_COUNTS.out.counts
    intron_chains = COMPARE_INTRON_CHAINS.out.results
    transcript_ends = COMPARE_TRANSCRIPT_ENDS.out.results
    utr_structure = COMPARE_UTR_STRUCTURE.out.results
    structural_summary = SUMMARISE_STRUCTURAL_COMPARISON.out.summary
    disagreement_loci = SUMMARISE_STRUCTURAL_COMPARISON.out.disagreements
    metric_status = METRIC_STATUS.out.status
    pair_report = ASSEMBLE_PAIR_REPORT.out.report
    evidence = ANNOTATE_MODEL_EVIDENCE.out.annotated
    gffcompare = gffcompare_summary
    versions = COLLECT_SOFTWARE_VERSIONS.out.versions
}
