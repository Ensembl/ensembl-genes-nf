include { TYPED_AXIS } from '../modules/typed_axis.nf'

workflow OBSERVED_EVIDENCE {
    take:
    peptides // [meta(peptide), peptide.json]

    main:
    requests = peptides.map { meta, file -> tuple(meta, 'observed_evidence', file) }
    TYPED_AXIS(requests)

    emit:
    axes = TYPED_AXIS.out.axis
    versions = TYPED_AXIS.out.versions
}
