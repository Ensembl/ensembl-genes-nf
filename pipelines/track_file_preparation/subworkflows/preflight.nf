include { VALIDATE_REFERENCE } from '../modules/validate_reference.nf'

workflow PREFLIGHT_TRACK_INPUTS {
    take:
    chrom_sizes
    assembly_release

    main:
    reference = VALIDATE_REFERENCE(chrom_sizes, assembly_release)

    emit:
    checked = reference.checked
    versions = reference.versions
}
