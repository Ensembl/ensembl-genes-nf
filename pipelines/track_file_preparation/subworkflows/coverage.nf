include { SAMTOOLS_QUICKCHECK } from '../modules/samtools_quickcheck.nf'
include { BAM_TO_BIGWIG } from '../modules/bam_to_bigwig.nf'
include { VALIDATE_BIGWIG } from '../modules/validate_bigwig.nf'
include { CHECKSUM_TRACK } from '../modules/checksum_track.nf'
include { WRITE_TRACK_RECORD } from '../modules/write_track_record.nf'

workflow PREPARE_COVERAGE_TRACKS {
    take:
    coverage_entities // tuple [meta, bam, bai, chrom_sizes]

    main:
    checked = SAMTOOLS_QUICKCHECK(coverage_entities)
    coverage = BAM_TO_BIGWIG(checked.checked)
    validated = VALIDATE_BIGWIG(coverage.bigwig)
    record_input = validated.validated.map { meta, track, source -> tuple(meta + [current_track: 'coverage'], track, source) }
    checksummed = CHECKSUM_TRACK(record_input)
    records = WRITE_TRACK_RECORD(checksummed.checksummed)

    emit:
    results = records.results
    versions = checked.versions.mix(coverage.versions).mix(records.versions)
}
