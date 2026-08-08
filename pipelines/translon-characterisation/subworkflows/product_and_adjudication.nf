include { PRODUCT_CLUSTERING } from '../modules/product_clustering.nf'
include { ADJUDICATE } from '../modules/adjudicate.nf'

workflow PRODUCT_AND_ADJUDICATION {
    take:
    characterised // [meta, characterised instances]

    main:
    // Both operations require the complete translon set; no earlier subworkflow collects.
    complete_set = characterised
        .map { meta, instance -> instance }
        .collect()
        .map { files -> tuple([id: 'translon-set'], files) }
    PRODUCT_CLUSTERING(complete_set)
    ADJUDICATE(complete_set, PRODUCT_CLUSTERING.out.clusters)

    emit:
    clusters = PRODUCT_CLUSTERING.out.clusters
    claims = ADJUDICATE.out.claims
    manifest = ADJUDICATE.out.manifest
    versions = PRODUCT_CLUSTERING.out.versions.concat(ADJUDICATE.out.versions)
}
