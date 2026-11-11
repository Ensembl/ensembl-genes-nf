nextflow.enable.dsl = 2

include { BUILD_CROSS_SPECIES } from '../modules/build_cross_species.nf'

workflow BUILD_CROSS_SPECIES_WORKFLOW {
    take:
    atlas
    orthology

    main:
    BUILD_CROSS_SPECIES(atlas.collect(), orthology, file("${projectDir}/bin/build_cross_species.py"))
}
