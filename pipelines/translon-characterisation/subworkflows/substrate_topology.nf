include { SUBSTRATE_TOPOLOGY } from '../modules/substrate_topology.nf'

workflow SUBSTRATE_TOPOLOGY_WORKFLOW {
    take:
    intervals  // [meta, intervals]
    annotation // [meta, gencode_gff3]
    genome // [meta, genome FASTA]

    main:
    SUBSTRATE_TOPOLOGY(intervals, annotation, genome)

    emit:
    instances = SUBSTRATE_TOPOLOGY.out.instances
    unhosted = SUBSTRATE_TOPOLOGY.out.unhosted
    versions = SUBSTRATE_TOPOLOGY.out.versions
}
