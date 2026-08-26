include { TYPED_AXIS } from '../modules/typed_axis.nf'
include { REGULATORY_GEOMETRY } from '../modules/regulatory_geometry.nf'

workflow INSTANCE_CONTEXT {
    take:
    instances // [meta, reconciled instances]

    main:
    REGULATORY_GEOMETRY(instances)
    requests = instances.combine(channel.of('expression_context'))
        .map { meta, file, axis -> tuple(meta, axis, file) }
    TYPED_AXIS(requests)

    emit:
    axes = TYPED_AXIS.out.axis.concat(REGULATORY_GEOMETRY.out.axis.map { meta, file -> tuple(meta, 'regulatory_geometry', file) })
    versions = TYPED_AXIS.out.versions.concat(REGULATORY_GEOMETRY.out.versions)
}
