include { GTF_TO_BIGBED } from '../modules/gtf_to_bigbed.nf'

workflow PREPARE_GENE_MODEL_TRACKS {
    take:
    gene_entities // tuple [meta, gtf, chrom_sizes]

    main:
    models = GTF_TO_BIGBED(gene_entities)

    emit:
    results = models.results
    versions = models.versions
}
