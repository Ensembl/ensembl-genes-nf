include { SAMTOOLS_QUICKCHECK } from '../modules/samtools_quickcheck.nf'
include { BAM_TO_BIGWIG } from '../modules/bam_to_bigwig.nf'
include { VALIDATE_BIGWIG } from '../../../modules/track/validate_bigwig.nf'
include { CHECKSUM_TRACK } from '../../../modules/track/checksum_track.nf'
include { WRITE_TRACK_RECORD } from '../modules/write_track_record.nf'

workflow PREPARE_COVERAGE_TRACKS {
    take:
    coverage_entities // tuple [meta, bam, bai, chrom_sizes]

    main:
    checked = SAMTOOLS_QUICKCHECK(coverage_entities)
    coverage = BAM_TO_BIGWIG(checked.checked)
    strand_tracks = coverage.bigwigs.flatMap { meta, forward, reverse, source ->
        [tuple(meta + [coverage_strand: 'forward', source_path: source.toString()], forward),
         tuple(meta + [coverage_strand: 'reverse', source_path: source.toString()], reverse)]
    }
    validated = VALIDATE_BIGWIG(strand_tracks)
    record_input = validated.validated.map { meta, track ->
        tuple(meta + [current_track: "coverage_${meta.coverage_strand}"], track)
    }
    checksummed = CHECKSUM_TRACK(record_input)
    records = WRITE_TRACK_RECORD(checksummed.checksummed)

    emit:
    results = records.results
    versions = checked.versions
        .mix(coverage.versions)
        .mix(validated.versions)
        .mix(checksummed.versions)
        .mix(records.versions)
}
