include { STAGE_TRANSCRIPT_INPUTS } from '../../modules/pairwise_annotation/stage_transcript_inputs.nf'
include { NORMALISE_TRANSCRIPT_PAIR } from '../../modules/pairwise_annotation/normalise_transcript_pair.nf'
include { COMPARE_TRANSCRIPT_FEATURES } from '../../modules/pairwise_annotation/compare_transcript_features.nf'

workflow PAIRWISE_TRANSCRIPT_COMPARISON {
    take:
    pairwise_samples_ch

    main:
    // Public input: tuple(meta, source_a, source_b, assembly_report,
    // evidence_a, evidence_b, reference_annotation).
    STAGE_TRANSCRIPT_INPUTS(pairwise_samples_ch, params.pairwise_cache_dir)
    staged = STAGE_TRANSCRIPT_INPUTS.out.staged

    NORMALISE_TRANSCRIPT_PAIR(staged.map { meta, source_a, source_b, _report, _evidence_a, _evidence_b, _reference ->
        tuple(meta, source_a, source_b)
    })
    models = NORMALISE_TRANSCRIPT_PAIR.out.models

    model_keyed = models.map { meta, source_a_models, source_b_models, audit ->
        tuple(meta.pair_id, meta, source_a_models, source_b_models, audit)
    }
    evidence_keyed = staged.map { meta, _source_a, _source_b, _report, evidence_a, evidence_b, _reference ->
        tuple(meta.pair_id, meta, evidence_a, evidence_b)
    }
    model_keyed.join(evidence_keyed).map { _pair_id, meta, source_a_models, source_b_models, _audit, _staged_meta, evidence_a, evidence_b ->
        tuple(meta, source_a_models, source_b_models, evidence_a, evidence_b)
    }.set { comparison_inputs }

    COMPARE_TRANSCRIPT_FEATURES(comparison_inputs)

    emit:
    staged_inputs = staged
    canonical_models = models.map { meta, source_a_models, source_b_models, _audit -> tuple(meta, source_a_models, source_b_models) }
    conversion_audits = models.map { meta, _source_a_models, _source_b_models, audit -> tuple(meta, audit) }
    model_presence = COMPARE_TRANSCRIPT_FEATURES.out.results.map { meta, model_presence, _feature_presence, _terminal_exons, _summary -> tuple(meta, model_presence) }
    feature_presence = COMPARE_TRANSCRIPT_FEATURES.out.results.map { meta, _model_presence, feature_presence, _terminal_exons, _summary -> tuple(meta, feature_presence) }
    terminal_exons = COMPARE_TRANSCRIPT_FEATURES.out.results.map { meta, _model_presence, _feature_presence, terminal_exons, _summary -> tuple(meta, terminal_exons) }
    structural_summary = COMPARE_TRANSCRIPT_FEATURES.out.results.map { meta, _model_presence, _feature_presence, _terminal_exons, summary -> tuple(meta, summary) }
}
