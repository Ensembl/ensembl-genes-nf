include { PREPARE_SPLICE_JUNCTIONS } from '../modules/splice_junction.nf'

workflow PREPARE_SPLICE_TRACKS {
    take:
    splice_entities // tuple [meta, sj_out_tab, bam, chrom_sizes]

    main:
    junctions = PREPARE_SPLICE_JUNCTIONS(splice_entities)

    emit:
    results = junctions.results
    versions = junctions.versions
}
