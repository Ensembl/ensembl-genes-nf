include { MERGE_AXES } from '../modules/merge_axes.nf'

workflow MERGE_TYPED_AXES {
    take:
    axes // [meta, axis_name, axis.jsonl]

    main:
    grouped = axes.map { meta, _axis, file -> tuple(meta, file) }.groupTuple()
    MERGE_AXES(grouped)

    emit:
    instances = MERGE_AXES.out.instances
    versions = MERGE_AXES.out.versions
}
