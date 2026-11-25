include { VALIDATE_REFERENCE } from '../modules/validate_reference.nf'
include { EXTRACT_BAM_HEADER } from '../modules/extract_bam_header.nf'
include { RESOLVE_BAM_REFERENCE } from '../modules/resolve_bam_reference.nf'

workflow PREFLIGHT_TRACK_INPUTS {
    take:
    chrom_sizes
    assembly_release
    bam_entities // tuple [meta, bam, chrom_sizes]

    main:
    reference = VALIDATE_REFERENCE(chrom_sizes, assembly_release,
        file("${projectDir}/bin/validate_reference.py", checkIfExists: true))
    bam_headers = EXTRACT_BAM_HEADER(bam_entities)
    resolved_bam_reference = RESOLVE_BAM_REFERENCE(bam_headers.header,
        file("${projectDir}/bin/resolve_bam_reference.py", checkIfExists: true))

    emit:
    checked = reference.checked
    bam_references = resolved_bam_reference.resolved
    versions = reference.versions.mix(bam_headers.versions).mix(resolved_bam_reference.versions)
}
