include { VALIDATE_PAIR_MANIFEST } from '../../modules/pairwise_annotation/validate_pair_manifest.nf'
include { STAGE_ANNOTATION_PAIR } from '../../modules/pairwise_annotation/stage_annotation_pair.nf'
include { NORMALISE_TRANSCRIPT_MODELS } from '../../modules/pairwise_annotation/normalise_transcript_models.nf'
include { AUDIT_ANNOTATION_CONVERSION } from '../../modules/pairwise_annotation/audit_annotation_conversion.nf'
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
    // Public boundary: tuple(meta, source_a, source_b, assembly_report,
    // evidence_a, evidence_b); every downstream channel retains meta.
    VALIDATE_PAIR_MANIFEST(pairwise_samples_ch)
    validated = VALIDATE_PAIR_MANIFEST.out.validated
    STAGE_ANNOTATION_PAIR(validated.map { meta, a, b, report, evidence_a, evidence_b, _validated -> tuple(meta, a, b, report, evidence_a, evidence_b) }, params.pairwise_cache_dir)
    staged = STAGE_ANNOTATION_PAIR.out.staged

    models_input = staged.map { meta, a, b, _report, _evidence_a, _evidence_b -> tuple(meta, a, b) }
    NORMALISE_TRANSCRIPT_MODELS(models_input)
    models = NORMALISE_TRANSCRIPT_MODELS.out.models
    AUDIT_ANNOTATION_CONVERSION(models)

    models_keyed = models.map { meta, source_a_models, source_b_models -> tuple(meta.pair_id, meta, source_a_models, source_b_models) }
    staged_keyed = staged.map { meta, _a, _b, _report, evidence_a, evidence_b -> tuple(meta.pair_id, meta, evidence_a, evidence_b) }
    models_keyed.join(staged_keyed).filter { _pair_id, _models_meta, _a_models, _b_models, _staged_meta, _evidence_a, _evidence_b -> _models_meta.evidence_a && _models_meta.evidence_b }.flatMap { _pair_id, meta, a_models, b_models, _staged_meta, evidence_a, evidence_b -> [tuple(meta, 'source_a', a_models, evidence_a), tuple(meta, 'source_b', b_models, evidence_b)] }.set { evidence_inputs }
    ANNOTATE_MODEL_EVIDENCE(evidence_inputs)

    COMPARE_INTRON_CHAINS(models)
    COMPARE_TRANSCRIPT_ENDS(models)
    COMPARE_UTR_STRUCTURE(models)

    intron_keyed = COMPARE_INTRON_CHAINS.out.results.map { meta, file -> tuple(meta.pair_id, meta, file) }
    ends_keyed = COMPARE_TRANSCRIPT_ENDS.out.results.map { meta, file -> tuple(meta.pair_id, meta, file) }
    utr_keyed = COMPARE_UTR_STRUCTURE.out.results.map { meta, file -> tuple(meta.pair_id, meta, file) }
    intron_keyed.join(ends_keyed).join(utr_keyed).map { _pair_id, meta, intron, _meta_ends, ends, _meta_utr, utr -> tuple(meta, intron, ends, utr) }.set { structural_summary_input }
    SUMMARISE_STRUCTURAL_COMPARISON(structural_summary_input)
    ASSEMBLE_PAIR_REPORT(SUMMARISE_STRUCTURAL_COMPARISON.out.summary.map { meta, summary -> tuple(meta, summary) })

    gffcompare_summary = channel.empty()
    gffcompare_versions = channel.empty()
    if (params.run_pairwise_gffcompare) {
        PREPARE_GFFCOMPARE_INPUTS(models)
        gffcompare_inputs = PREPARE_GFFCOMPARE_INPUTS.out.gtfs.flatMap { meta, source_a_gtf, source_b_gtf -> [
            tuple(meta, 'source_a_to_source_b', source_a_gtf, source_b_gtf),
            tuple(meta, 'source_b_to_source_a', source_b_gtf, source_a_gtf)
        ] }
        GFFCOMPARE(gffcompare_inputs)
        PARSE_GFFCOMPARE_TMAP(GFFCOMPARE.out.tmap)
        gffcompare_classes = PARSE_GFFCOMPARE_TMAP.out.classes.map { meta, _direction, classes -> tuple(meta.pair_id, meta, classes) }.groupTuple().map { _pair_id, metas, classes -> tuple(metas[0], classes) }
        SUMMARISE_GFFCOMPARE_MODELS(gffcompare_classes)
        gffcompare_summary = SUMMARISE_GFFCOMPARE_MODELS.out.summary
        gffcompare_versions = PREPARE_GFFCOMPARE_INPUTS.out.versions
            .mix(GFFCOMPARE.out.versions)
            .mix(PARSE_GFFCOMPARE_TMAP.out.versions)
            .mix(SUMMARISE_GFFCOMPARE_MODELS.out.versions)
    }

    status_input = pairwise_samples_ch.map { meta, _a, _b, _report, _ea, _eb ->
        def bed12 = meta.format_a == 'bed12' || meta.format_b == 'bed12'
        def entries = [
            ['intron_chain', 'enabled', 'canonical model comparison'],
            ['transcript_ends', 'enabled', 'canonical model comparison'],
            ['utr_structure', 'enabled', 'available fields are reported per model'],
            ['gene_presence', bed12 ? 'skipped' : 'enabled', bed12 ? 'BED12 has no biological gene contract' : 'annotation gene records'],
            ['coding_integrity', bed12 ? 'skipped' : 'enabled', bed12 ? 'BED12 has no CDS contract' : 'annotation CDS records'],
            ['gff_feature_metrics', bed12 ? 'skipped' : 'enabled', bed12 ? 'BED12 has no gene/CDS feature contract' : 'annotation feature records'],
            ['gffcompare', params.run_pairwise_gffcompare ? 'enabled' : 'skipped', params.run_pairwise_gffcompare ? 'requested' : 'disabled by parameter'],
        ]
        tuple(meta, entries)
    }
    METRIC_STATUS(status_input)

    pairwise_versions = VALIDATE_PAIR_MANIFEST.out.versions
        .mix(STAGE_ANNOTATION_PAIR.out.versions)
        .mix(NORMALISE_TRANSCRIPT_MODELS.out.versions)
        .mix(AUDIT_ANNOTATION_CONVERSION.out.versions)
        .mix(COMPARE_INTRON_CHAINS.out.versions)
        .mix(COMPARE_TRANSCRIPT_ENDS.out.versions)
        .mix(COMPARE_UTR_STRUCTURE.out.versions)
        .mix(SUMMARISE_STRUCTURAL_COMPARISON.out.versions)
        .mix(METRIC_STATUS.out.versions)
        .mix(ASSEMBLE_PAIR_REPORT.out.versions)
        .mix(ANNOTATE_MODEL_EVIDENCE.out.versions)
        .mix(gffcompare_versions)
        .collect()
    COLLECT_SOFTWARE_VERSIONS(pairwise_versions)

    // Preserve established GFF3 metrics for non-BED12 pairs. Structural
    // metrics above remain the only metrics enabled for BED12 inputs.
    legacy = staged.filter { meta, _a, _b, _report, _ea, _eb -> meta.format_a != 'bed12' && meta.format_b != 'bed12' }
    if (params.run_pairwise_legacy_metrics) {
        PAIR_GENES(legacy.map { meta, a, b, report, _ea, _eb -> tuple(meta, a, b, report) })
        rbh_inputs = legacy.map { meta, a, b, _report, _ea, _eb -> tuple(meta, a, b) }.join(PAIR_GENES.out.rbh)
        if (params.run_pairwise_transcript_concordance) {
            PAIRWISE_TRANSCRIPT_CONCORDANCE(rbh_inputs)
        }
        if (params.run_pairwise_coding_integrity) {
            PAIRWISE_CODING_INTEGRITY(rbh_inputs)
        }
        if (params.run_pairwise_gene_presence) {
            PAIRWISE_GENE_PRESENCE(legacy.map { meta, a, b, _report, _ea, _eb -> tuple(meta, a, b) }, '')
        }
        if (params.run_pairwise_multi_mapping) {
            PAIRWISE_MULTI_MAPPING(PAIR_GENES.out.all_pairs)
        }
        if (params.run_pairwise_gff_feature_metrics) {
            PAIRWISE_GFF_FEATURE_METRICS(legacy.map { meta, a, b, _report, _ea, _eb -> tuple(meta, a, b) })
        }
        if (params.run_pairwise_transcript_counts) {
            PAIRWISE_COUNT_GFF_TRANSCRIPTS(legacy.flatMap { meta, a, b, _report, _ea, _eb -> [tuple(meta, 'source_a', a), tuple(meta, 'source_b', b)] })
        }
    }

    emit:
    staged_inputs = staged
    models = models
    conversion_audits = AUDIT_ANNOTATION_CONVERSION.out.audit
    id_maps = AUDIT_ANNOTATION_CONVERSION.out.id_map
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
