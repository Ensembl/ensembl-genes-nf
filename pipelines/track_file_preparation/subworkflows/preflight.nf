include { VALIDATE_REFERENCE } from '../modules/validate_reference.nf'
include { RESOLVE_BAM_REFERENCE } from '../modules/resolve_bam_reference.nf'

workflow PREFLIGHT_TRACK_INPUTS {
    take:
    chrom_sizes
    assembly_release
    bam_entities // tuple [meta, bam, chrom_sizes]

    main:
    reference = VALIDATE_REFERENCE(chrom_sizes, assembly_release)
    resolved_bam_reference = RESOLVE_BAM_REFERENCE(bam_entities)

    emit:
    checked = reference.checked
    bam_references = resolved_bam_reference.resolved
    versions = reference.versions.mix(resolved_bam_reference.versions)
}
