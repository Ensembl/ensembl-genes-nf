include { VALIDATE_STAR_JUNCTIONS } from '../modules/validate_star_junctions.nf'
include { STAR_JUNCTIONS_TO_BED } from '../modules/sj_to_bed.nf'
include { SORT_JUNCTION_BED } from '../modules/sort_junction_bed.nf'
include { JUNCTION_BED_TO_BIGBED } from '../../../modules/track/junction_bed_to_bigbed.nf'
include { VALIDATE_SPLICE_BIGBED } from '../../../modules/track/validate_splice_bigbed.nf'
include { CHECKSUM_TRACK } from '../../../modules/track/checksum_track.nf'
include { WRITE_TRACK_RECORD } from '../modules/write_track_record.nf'

workflow PREPARE_SPLICE_TRACKS {
    take:
    splice_entities // tuple [meta, sj_out_tab, bam, chrom_sizes]

    main:
    splice_inputs = splice_entities.map { meta, sj, bam, chrom_sizes ->
        tuple(meta + [source_path: sj.toString()], sj, bam, chrom_sizes)
    }
    supplied_validated = VALIDATE_STAR_JUNCTIONS(splice_inputs)
    bed = STAR_JUNCTIONS_TO_BED(supplied_validated.validated, file("${projectDir}/bin/sj_to_bed.py", checkIfExists: true))
    sorted = SORT_JUNCTION_BED(bed.bed)
    bigbed_input = sorted.sorted.map { meta, junction_bed, _source, chrom_sizes -> tuple(meta, junction_bed, chrom_sizes) }
    bigbed = JUNCTION_BED_TO_BIGBED(bigbed_input, file("${projectDir}/resources/splice_junction.as", checkIfExists: true))
    validated = VALIDATE_SPLICE_BIGBED(bigbed.bigbed)
    record_input = validated.validated.map { meta, track -> tuple(meta + [current_track: 'splice_junction'], track) }
    checksummed = CHECKSUM_TRACK(record_input)
    records = WRITE_TRACK_RECORD(checksummed.checksummed)

    emit:
    results = records.results
    versions = bigbed.versions
        .mix(validated.versions)
        .mix(checksummed.versions)
        .mix(records.versions)
}
