include { TYPED_AXIS } from '../modules/typed_axis.nf'
include { PEPTIDE_UNIQUENESS } from '../modules/peptide_uniqueness.nf'
include { FANBACK_PEPTIDE_AXIS } from '../modules/fanback_peptide_axis.nf'

workflow PEPTIDE_AXES {
    take:
    peptides // [meta(peptide), peptide.json]
    proteome // [meta, GENCODE proteome FASTA]
    instances // [meta(instance), reconciled instances]

    main:
    PEPTIDE_UNIQUENESS(peptides, proteome)
    FANBACK_PEPTIDE_AXIS(instances, PEPTIDE_UNIQUENESS.out.results)
    requests = instances.combine(channel.of('peptide_features', 'structure', 'detectability'))
        .map { meta, file, axis -> tuple(meta, axis, file) }
    TYPED_AXIS(requests)

    emit:
    axes = TYPED_AXIS.out.axis.concat(FANBACK_PEPTIDE_AXIS.out.axis.map { meta, file -> tuple(meta, 'peptide_uniqueness', file) })
    versions = TYPED_AXIS.out.versions.concat(PEPTIDE_UNIQUENESS.out.versions).concat(FANBACK_PEPTIDE_AXIS.out.versions)
}
