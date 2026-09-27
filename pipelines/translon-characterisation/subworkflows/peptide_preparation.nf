include { PEPTIDE_PRODUCTS } from '../modules/peptide_products.nf'

workflow PEPTIDE_PREPARATION {
    take:
    instances // [meta, reconciled instances]

    main:
    PEPTIDE_PRODUCTS(instances)

    emit:
    peptides = PEPTIDE_PRODUCTS.out.peptides
    fanback = PEPTIDE_PRODUCTS.out.fanback
    versions = PEPTIDE_PRODUCTS.out.versions
}
