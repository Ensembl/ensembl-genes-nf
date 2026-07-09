/*
 * PANGENOME_PROJECTION
 * Project reference gene annotation onto one or more target assemblies.
 *
 * A swappable whole-genome aligner (MINIMAP2_WGA) produces a PAF that the hpp
 * projection stages consume. The hpp stages run in the monolith's real order:
 *   synteny -> project -> validate -> rescue -> refine -> finalise
 * Each target assembly flows through the DAG independently (per-target parallelism).
 */

include { MINIMAP2_WGA       } from '../modules/minimap2_wga.nf'
include { BUILD_SYNTENY       } from '../modules/build_synteny.nf'
include { PROJECT_FEATURES    } from '../modules/project_features.nf'
include { VALIDATE_MODELS     } from '../modules/validate_models.nf'
include { RESCUE_PROJECTIONS  } from '../modules/rescue_projections.nf'
include { REFINE_MODELS       } from '../modules/refine_models.nf'
include { FINALISE            } from '../modules/finalise.nf'

workflow PANGENOME_PROJECTION {
    take:
    targets        // channel: [ val(meta), path(target_fasta) ]  (meta.id = target assembly id)
    reference      // path: reference genome FASTA (shared)
    annotation     // path: reference annotation GFF3 (shared)

    main:
    ch_versions = Channel.empty()

    // 1. Whole-genome alignment (swappable front-end) -> PAF
    MINIMAP2_WGA(targets, reference)
    ch_versions = ch_versions.mix(MINIMAP2_WGA.out.versions.first())

    // 2. Parse synteny from the PAF (+ gap fill, sex chroms).
    //    Join the PAF back with its target assembly: [ meta, paf, target_fasta ]
    ch_build_in = MINIMAP2_WGA.out.paf.join(targets)
    BUILD_SYNTENY(ch_build_in, reference, annotation)
    ch_versions = ch_versions.mix(BUILD_SYNTENY.out.versions.first())

    // 3-6. Sequential hpp stages threading the state directory. Each FASTA-consuming
    //      stage gets the target assembly staged in (joined by meta) plus the shared
    //      reference, so FastaHandler resolves them in each separate workdir.
    PROJECT_FEATURES(BUILD_SYNTENY.out.state.join(targets), reference)
    VALIDATE_MODELS(PROJECT_FEATURES.out.state.join(targets), reference)
    RESCUE_PROJECTIONS(VALIDATE_MODELS.out.state.join(targets), reference)
    REFINE_MODELS(RESCUE_PROJECTIONS.out.state.join(targets), reference)
    ch_versions = ch_versions
        .mix(PROJECT_FEATURES.out.versions.first())
        .mix(VALIDATE_MODELS.out.versions.first())
        .mix(RESCUE_PROJECTIONS.out.versions.first())
        .mix(REFINE_MODELS.out.versions.first())

    // 7. Finalise -> GFF3 + stats
    FINALISE(REFINE_MODELS.out.state)
    ch_versions = ch_versions.mix(FINALISE.out.versions.first())

    emit:
    gff      = FINALISE.out.gff       // channel: [ meta, path(mapped.gff3) ]
    stats    = FINALISE.out.stats     // channel: [ meta, path(stats.json) ]
    versions = ch_versions            // channel: path(versions.yml)
}
