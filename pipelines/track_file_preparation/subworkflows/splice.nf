include { VALIDATE_STAR_JUNCTIONS } from '../modules/validate_star_junctions.nf'
include { STAR_JUNCTIONS_TO_BED } from '../modules/sj_to_bed.nf'
include { SORT_JUNCTION_BED } from '../modules/sort_junction_bed.nf'
include { JUNCTION_BED_TO_BIGBED } from '../modules/junction_bed_to_bigbed.nf'
include { VALIDATE_SPLICE_BIGBED } from '../modules/validate_splice_bigbed.nf'
include { CHECKSUM_TRACK } from '../modules/checksum_track.nf'
include { WRITE_TRACK_RECORD } from '../modules/write_track_record.nf'

workflow PREPARE_SPLICE_TRACKS {
    take:
    splice_entities // tuple [meta, sj_out_tab, bam, chrom_sizes]

    main:
    supplied_validated = VALIDATE_STAR_JUNCTIONS(splice_entities)
    bed = STAR_JUNCTIONS_TO_BED(supplied_validated.validated, file("${projectDir}/bin/sj_to_bed.py", checkIfExists: true))
    sorted = SORT_JUNCTION_BED(bed.bed)
    bigbed = JUNCTION_BED_TO_BIGBED(sorted.sorted, file("${projectDir}/resources/splice_junction.as", checkIfExists: true))
    validated = VALIDATE_SPLICE_BIGBED(bigbed.bigbed)
    record_input = validated.validated.map { meta, track, source -> tuple(meta + [current_track: 'splice_junction'], track, source) }
    checksummed = CHECKSUM_TRACK(record_input)
    records = WRITE_TRACK_RECORD(checksummed.checksummed)

    emit:
    results = records.results
    versions = records.versions
}
