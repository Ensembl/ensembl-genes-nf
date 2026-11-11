nextflow.enable.dsl = 2

include { BUILD_SPECIES_ATLAS } from '../modules/build_species_atlas.nf'

workflow BUILD_SPECIES_ATLAS_WORKFLOW {
    take:
    ctss
    annotations

    main:
    ch_grouped = ctss.map { meta, bed, bedgraph -> tuple(meta.assembly, meta, bed, bedgraph) }
        .groupTuple()
        .map { _assembly, metas, beds, bedgraphs -> tuple(_assembly, metas[0], metas, beds, bedgraphs) }
    ch_atlas_input = ch_grouped.join(annotations)
        .map { _assembly, meta, metas, beds, bedgraphs, gtf -> tuple(meta, metas, beds, bedgraphs, gtf) }
    ch_atlas = BUILD_SPECIES_ATLAS(ch_atlas_input, file("${projectDir}/bin/build_species_atlas.R")).species_atlas

    emit:
    species_atlas = ch_atlas
}
