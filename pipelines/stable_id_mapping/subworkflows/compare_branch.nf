include { STAGE_SPECIES_INPUTS } from '../modules/stage_species_inputs.nf'
include { LIFTON_PROJECTION } from '../modules/lifton_projection.nf'
include { STRUCTURAL_MATCHING } from '../modules/structural_matching.nf'
include { COMPARE_ANNOTATIONS } from '../modules/compare_annotations.nf'


workflow COMPARE_BRANCH {
    take:
    comparison_ch

    main:
    rules_config_ch = Channel.value(
        file(params.rules_config, checkIfExists: true)
    )

    staged_comparison_ch = comparison_ch.map {
        comparison_name, ref_fasta, ref_gff, target_fasta, target_gff ->
        tuple(
            comparison_name,
            ref_fasta,
            ref_gff,
            target_fasta,
            target_gff,
            0,
            'COMPAREG:1-999999999999',
            'COMPARET:1-999999999999',
            'COMPAREP:1-999999999999',
        )
    }

    STAGE_SPECIES_INPUTS(staged_comparison_ch)
    LIFTON_PROJECTION(STAGE_SPECIES_INPUTS.out.staged_inputs)
    STRUCTURAL_MATCHING(LIFTON_PROJECTION.out.projected, rules_config_ch)
    COMPARE_ANNOTATIONS(STRUCTURAL_MATCHING.out.matches, rules_config_ch)

    emit:
    relationships = COMPARE_ANNOTATIONS.out.relationships
}
