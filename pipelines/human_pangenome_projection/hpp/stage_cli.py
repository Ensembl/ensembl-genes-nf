"""Per-stage CLI commands for the Nextflow-orchestrated pipeline.

Each command hydrates a :class:`~hpp.pipeline.PangenomeMappingPipeline` from a
serialized state directory, runs the *unmodified* monolith stage method, and
writes the updated state back out. Because the stage logic is the monolith's own
methods (never reimplemented), running the stages in sequence is byte-identical
to ``hpp map``.

Stage order mirrors the monolith's real control flow:
    build-synteny -> project-features -> validate-models
    -> rescue-projections -> refine-models -> finalise
"""

import logging
import sys
from pathlib import Path

import click

from hpp.config import MappingConfig
from hpp.pipeline import PangenomeMappingPipeline
from hpp.state import dump_state, hydrate_pipeline

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Shared configuration options (identical set used by `map` and `build-synteny`)
# ---------------------------------------------------------------------------

_CONFIG_OPTIONS = [
    click.option("--ref-fasta", "-r", required=True, type=click.Path(exists=True),
                 help="Reference genome FASTA file"),
    click.option("--ref-gff", "-g", required=True, type=click.Path(exists=True),
                 help="Reference annotation GFF3 file"),
    click.option("--target-fasta", "-t", required=True, type=click.Path(exists=True),
                 help="Target genome FASTA file"),
    click.option("--output-gff", "-o", required=True, type=click.Path(),
                 help="Output GFF3 file for mapped annotations"),
    click.option("--output-stats", "-s", type=click.Path(),
                 help="Output JSON file for mapping statistics"),
    click.option("--chromosomes", "-c", multiple=True,
                 help="Limit to specific chromosomes (repeatable)"),
    click.option("--threads", "-j", default=4, type=int,
                 help="Number of threads for alignment (default: 4)"),
    click.option("--preset-profile", default="none",
                 type=click.Choice(["none", "ultra-close", "less-close"]),
                 help="Production parameter profile (default: none)"),
    click.option("--min-identity", default=0.95, type=float),
    click.option("--min-block-length", default=10000, type=int),
    click.option("--min-mapq", default=10, type=int),
    click.option("--cnv-max-total-copies", default=3, type=int),
    click.option("--cnv-min-expansion-coverage", default=0.60, type=float),
    click.option("--cnv-min-score-ratio", default=0.85, type=float),
    click.option("--cnv-max-locus-overlap", default=0.30, type=float),
    click.option("--cnv-group-min-reciprocal-overlap", default=0.50, type=float),
    click.option("--cnv-ambiguity-distance-bp", default=10000, type=int),
    click.option("--cnv-ambiguity-score-delta", default=0.01, type=float),
    click.option("--disable-paralog-reassignment", is_flag=True),
    click.option("--paralog-overlap-threshold", default=0.50, type=float),
    click.option("--paralog-max-locus-candidates", default=8, type=int),
    click.option("--paralog-max-rounds", default=2, type=int),
    click.option("--disable-missed-locus-recovery", is_flag=True),
    click.option("--missed-locus-max-genes", default=1000, type=int),
    click.option("--missed-locus-search-padding-bp", default=200000, type=int),
    click.option("--missed-locus-max-window-bp", default=800000, type=int),
    click.option("--missed-locus-min-identity", default=0.85, type=float),
    click.option("--missed-locus-min-coverage", default=0.60, type=float),
    click.option("--disable-boundary-refinement", is_flag=True),
    click.option("--boundary-refine-min-coverage", default=0.98, type=float),
    click.option("--boundary-refine-window-bp", default=500, type=int),
    click.option("--boundary-refine-anchor-bp", default=80, type=int),
    click.option("--boundary-refine-max-unaligned-bp", default=8, type=int),
    click.option("--disable-validation-rescue", is_flag=True),
    click.option("--validation-rescue-window-bp", default=1200, type=int),
    click.option("--validation-rescue-anchor-bp", default=120, type=int),
    click.option("--validation-rescue-max-genes", default=2000, type=int),
    click.option("--disable-internal-realign-rescue", is_flag=True),
    click.option("--internal-realign-max-genes", default=3000, type=int),
    click.option("--internal-realign-ref-flank-bp", default=3000, type=int),
    click.option("--internal-realign-target-flank-bp", default=10000, type=int),
    click.option("--internal-realign-min-path-coverage", default=0.98, type=float),
    click.option("--internal-realign-max-internal-gap-bp", default=6, type=int),
    click.option("--internal-realign-trigger-indel-jump-bp", default=6, type=int),
    click.option("--internal-realign-fallback-backend", default="edlib",
                 type=click.Choice(["edlib"])),
    click.option("--disable-internal-realign-accept-only-if-improved", is_flag=True),
    click.option("--disable-splice-shift-rescue", is_flag=True),
    click.option("--splice-shift-max-bp", default=6, type=int),
    click.option("--disable-codon-frame-rescue", is_flag=True),
    click.option("--codon-rescue-max-bp", default=6, type=int),
    click.option("--disable-protein-qc", is_flag=True),
    click.option("--protein-qc-min-identity", default=0.90, type=float),
    click.option("--protein-qc-min-coverage", default=0.90, type=float),
    click.option("--allow-secondary/--primary-only", default=False),
    click.option("--keep-temp", is_flag=True),
    click.option("--verbose", "-v", is_flag=True),
]


def config_options(func):
    """Apply the full shared option set to a click command."""
    for option in reversed(_CONFIG_OPTIONS):
        func = option(func)
    return func


def resolve_config(p: dict) -> MappingConfig:
    """Build a :class:`MappingConfig` from click params (incl. preset profiles).

    This is the exact resolution the monolith ``map`` command performs, factored
    out so every entry point produces an identical config.
    """
    chr_set = set(p["chromosomes"]) if p.get("chromosomes") else None

    min_identity = p["min_identity"]
    cnv_min_expansion_coverage = p["cnv_min_expansion_coverage"]
    cnv_min_score_ratio = p["cnv_min_score_ratio"]
    paralog_overlap_threshold = p["paralog_overlap_threshold"]
    missed_locus_min_identity = p["missed_locus_min_identity"]
    missed_locus_min_coverage = p["missed_locus_min_coverage"]
    boundary_refine_window_bp = p["boundary_refine_window_bp"]
    validation_rescue_window_bp = p["validation_rescue_window_bp"]
    missed_locus_search_padding_bp = p["missed_locus_search_padding_bp"]

    profile = p.get("preset_profile", "none")
    if profile == "ultra-close":
        min_identity = max(min_identity, 0.98)
        cnv_min_expansion_coverage = max(cnv_min_expansion_coverage, 0.70)
        cnv_min_score_ratio = max(cnv_min_score_ratio, 0.90)
        paralog_overlap_threshold = max(paralog_overlap_threshold, 0.55)
        missed_locus_min_identity = max(missed_locus_min_identity, 0.90)
        missed_locus_min_coverage = max(missed_locus_min_coverage, 0.70)
    elif profile == "less-close":
        min_identity = min(min_identity, 0.93)
        cnv_min_expansion_coverage = min(cnv_min_expansion_coverage, 0.55)
        cnv_min_score_ratio = min(cnv_min_score_ratio, 0.80)
        boundary_refine_window_bp = max(boundary_refine_window_bp, 800)
        validation_rescue_window_bp = max(validation_rescue_window_bp, 1600)
        missed_locus_search_padding_bp = max(missed_locus_search_padding_bp, 300000)

    return MappingConfig(
        ref_fasta=Path(p["ref_fasta"]),
        ref_gff=Path(p["ref_gff"]),
        target_fasta=Path(p["target_fasta"]),
        output_gff=Path(p["output_gff"]),
        output_stats=Path(p["output_stats"]) if p.get("output_stats") else None,
        chromosomes=chr_set,
        threads=p["threads"],
        min_identity=min_identity,
        min_block_length=p["min_block_length"],
        min_mapq=p["min_mapq"],
        primary_only=not p["allow_secondary"],
        cnv_max_total_copies=p["cnv_max_total_copies"],
        cnv_min_expansion_coverage=cnv_min_expansion_coverage,
        cnv_min_score_ratio=cnv_min_score_ratio,
        cnv_max_locus_overlap=p["cnv_max_locus_overlap"],
        cnv_group_min_reciprocal_overlap=p["cnv_group_min_reciprocal_overlap"],
        cnv_ambiguity_distance_bp=p["cnv_ambiguity_distance_bp"],
        cnv_ambiguity_score_delta=p["cnv_ambiguity_score_delta"],
        enable_paralog_reassignment=not p["disable_paralog_reassignment"],
        paralog_overlap_threshold=paralog_overlap_threshold,
        paralog_max_locus_candidates=p["paralog_max_locus_candidates"],
        paralog_max_rounds=p["paralog_max_rounds"],
        enable_missed_locus_recovery=not p["disable_missed_locus_recovery"],
        missed_locus_max_genes=p["missed_locus_max_genes"],
        missed_locus_search_padding_bp=missed_locus_search_padding_bp,
        missed_locus_max_window_bp=p["missed_locus_max_window_bp"],
        missed_locus_min_identity=missed_locus_min_identity,
        missed_locus_min_coverage=missed_locus_min_coverage,
        enable_boundary_refinement=not p["disable_boundary_refinement"],
        boundary_refine_min_coverage=p["boundary_refine_min_coverage"],
        boundary_refine_window_bp=boundary_refine_window_bp,
        boundary_refine_anchor_bp=p["boundary_refine_anchor_bp"],
        boundary_refine_max_unaligned_bp=p["boundary_refine_max_unaligned_bp"],
        enable_validation_rescue=not p["disable_validation_rescue"],
        validation_rescue_window_bp=validation_rescue_window_bp,
        validation_rescue_anchor_bp=p["validation_rescue_anchor_bp"],
        validation_rescue_max_genes=p["validation_rescue_max_genes"],
        enable_internal_realign_rescue=not p["disable_internal_realign_rescue"],
        internal_realign_max_genes=p["internal_realign_max_genes"],
        internal_realign_ref_flank_bp=p["internal_realign_ref_flank_bp"],
        internal_realign_target_flank_bp=p["internal_realign_target_flank_bp"],
        internal_realign_min_path_coverage=p["internal_realign_min_path_coverage"],
        internal_realign_max_internal_gap_bp=p["internal_realign_max_internal_gap_bp"],
        internal_realign_trigger_indel_jump_bp=p["internal_realign_trigger_indel_jump_bp"],
        internal_realign_fallback_backend=p["internal_realign_fallback_backend"],
        internal_realign_accept_only_if_improved=(
            not p["disable_internal_realign_accept_only_if_improved"]
        ),
        enable_splice_shift_rescue=not p["disable_splice_shift_rescue"],
        splice_shift_max_bp=p["splice_shift_max_bp"],
        enable_codon_frame_rescue=not p["disable_codon_frame_rescue"],
        codon_rescue_max_bp=p["codon_rescue_max_bp"],
        enable_protein_qc=not p["disable_protein_qc"],
        protein_qc_min_identity=p["protein_qc_min_identity"],
        protein_qc_min_coverage=p["protein_qc_min_coverage"],
        keep_temp=p["keep_temp"],
    )


def fasta_options(func):
    """Add --ref-fasta/--target-fasta overrides (staged copies in a process workdir)."""
    func = click.option("--ref-fasta", "ref_fasta", type=click.Path(exists=True), default=None,
                        help="Reference FASTA (overrides the path stored in run config)")(func)
    func = click.option("--target-fasta", "target_fasta", type=click.Path(exists=True), default=None,
                        help="Target FASTA (overrides the path stored in run config)")(func)
    return func


def _apply_fasta_overrides(pipeline, ref_fasta, target_fasta) -> None:
    """Point the run config at locally-staged FASTAs when provided.

    Staged stages run in separate workdirs, so the absolute/basename paths captured
    at build-synteny time may not resolve. The Nextflow modules stage the FASTAs in
    and pass them here so FastaHandler opens the right files.
    """
    if ref_fasta:
        pipeline.config.ref_fasta = Path(ref_fasta)
    if target_fasta:
        pipeline.config.target_fasta = Path(target_fasta)


def _setup_logging(verbose: bool) -> None:
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)


def _fail(exc: Exception, verbose: bool) -> None:
    logger.error(f"Stage failed: {exc}")
    if verbose:
        import traceback
        traceback.print_exc()
    sys.exit(1)


# ---------------------------------------------------------------------------
# Stage commands
# ---------------------------------------------------------------------------

@click.command("build-synteny")
@config_options
@click.option("--paf", type=click.Path(exists=True), default=None,
              help="Pre-computed cs-tagged PAF from a standalone WGA module. "
                   "When given, the main whole-genome alignment is not run here.")
@click.option("--out-state", required=True, type=click.Path(),
              help="Output state directory for downstream stages")
def build_synteny(**params):
    """Stage 1: load inputs, detect synteny, detect sex chromosomes."""
    _setup_logging(params["verbose"])
    try:
        config = resolve_config(params)
        pipeline = PangenomeMappingPipeline(config)
        if params.get("paf"):
            pipeline.external_paf = Path(params["paf"])
        pipeline._load_inputs()
        pipeline._run_synteny_detection()
        pipeline._detect_sex_chromosomes()
        dump_state(pipeline, params["out_state"])
        click.echo(f"build-synteny: wrote state to {params['out_state']}")
    except Exception as exc:  # noqa: BLE001
        _fail(exc, params["verbose"])


@click.command("project-features")
@click.option("--in-state", required=True, type=click.Path(exists=True))
@click.option("--out-state", required=True, type=click.Path())
@fasta_options
@click.option("--verbose", "-v", is_flag=True)
def project_features(in_state, out_state, ref_fasta, target_fasta, verbose):
    """Stage 2: project gene features through synteny (incl. paralog + recovery)."""
    _setup_logging(verbose)
    try:
        pipeline = hydrate_pipeline(in_state)
        _apply_fasta_overrides(pipeline, ref_fasta, target_fasta)
        pipeline._run_feature_mapping()
        dump_state(pipeline, out_state)
        click.echo(f"project-features: wrote state to {out_state}")
    except Exception as exc:  # noqa: BLE001
        _fail(exc, verbose)


@click.command("validate-models")
@click.option("--in-state", required=True, type=click.Path(exists=True))
@click.option("--out-state", required=True, type=click.Path())
@fasta_options
@click.option("--verbose", "-v", is_flag=True)
def validate_models(in_state, out_state, ref_fasta, target_fasta, verbose):
    """Stage 3: initial structural validation of projected models."""
    _setup_logging(verbose)
    try:
        pipeline = hydrate_pipeline(in_state)
        _apply_fasta_overrides(pipeline, ref_fasta, target_fasta)
        pipeline._validate_initial()
        dump_state(pipeline, out_state)
        click.echo(f"validate-models: wrote state to {out_state}")
    except Exception as exc:  # noqa: BLE001
        _fail(exc, verbose)


@click.command("rescue-projections")
@click.option("--in-state", required=True, type=click.Path(exists=True))
@click.option("--out-state", required=True, type=click.Path())
@fasta_options
@click.option("--verbose", "-v", is_flag=True)
def rescue_projections(in_state, out_state, ref_fasta, target_fasta, verbose):
    """Stage 4: rescue passes, re-validation, and protein QC."""
    _setup_logging(verbose)
    try:
        pipeline = hydrate_pipeline(in_state)
        _apply_fasta_overrides(pipeline, ref_fasta, target_fasta)
        pipeline._run_rescues_and_protein()
        dump_state(pipeline, out_state)
        click.echo(f"rescue-projections: wrote state to {out_state}")
    except Exception as exc:  # noqa: BLE001
        _fail(exc, verbose)


@click.command("refine-models")
@click.option("--in-state", required=True, type=click.Path(exists=True))
@click.option("--out-state", required=True, type=click.Path())
@fasta_options
@click.option("--verbose", "-v", is_flag=True)
def refine_models(in_state, out_state, ref_fasta, target_fasta, verbose):
    """Stage 5: conflict resolution / refinement and audit traces."""
    _setup_logging(verbose)
    try:
        pipeline = hydrate_pipeline(in_state)
        _apply_fasta_overrides(pipeline, ref_fasta, target_fasta)
        pipeline._run_refinement()
        dump_state(pipeline, out_state)
        click.echo(f"refine-models: wrote state to {out_state}")
    except Exception as exc:  # noqa: BLE001
        _fail(exc, verbose)


@click.command("finalise")
@click.option("--in-state", required=True, type=click.Path(exists=True))
@click.option("--output-gff", "-o", type=click.Path(),
              help="Override output GFF3 path (defaults to run config)")
@click.option("--output-stats", "-s", type=click.Path(),
              help="Override output stats JSON path (defaults to run config)")
@click.option("--verbose", "-v", is_flag=True)
def finalise(in_state, output_gff, output_stats, verbose):
    """Stage 6: statistics, synteny analysis, and final output generation."""
    _setup_logging(verbose)
    try:
        pipeline = hydrate_pipeline(in_state)
        if output_gff:
            pipeline.config.output_gff = Path(output_gff)
        if output_stats:
            pipeline.config.output_stats = Path(output_stats)
        pipeline._calculate_statistics()
        pipeline._run_synteny_analysis()
        # Timing fields are run-specific and stripped by the equivalence harness.
        pipeline.stats.total_time_seconds = 0.0
        pipeline._generate_outputs()
        click.echo(f"finalise: wrote {pipeline.config.output_gff}")
    except Exception as exc:  # noqa: BLE001
        _fail(exc, verbose)


STAGE_COMMANDS = [
    build_synteny,
    project_features,
    validate_models,
    rescue_projections,
    refine_models,
    finalise,
]
