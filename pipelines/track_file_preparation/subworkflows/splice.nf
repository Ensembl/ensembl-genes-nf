include { VALIDATE_STAR_JUNCTIONS } from '../modules/validate_star_junctions.nf'
include { EXTRACT_JUNCTIONS } from '../modules/extract_junctions.nf'
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
    supplied = splice_entities.filter { meta, sj, bam, chrom_sizes -> !meta.derive_junctions }
    derived = splice_entities.filter { meta, sj, bam, chrom_sizes -> meta.derive_junctions }
    supplied_validated = VALIDATE_STAR_JUNCTIONS(supplied)
    derived_junctions = EXTRACT_JUNCTIONS(derived.map { meta, sj, bam, chrom_sizes -> tuple(meta, bam, chrom_sizes) })
    junctions = supplied_validated.validated.mix(derived_junctions.junctions)
    bed = STAR_JUNCTIONS_TO_BED(junctions)
    sorted = SORT_JUNCTION_BED(bed.bed)
    bigbed = JUNCTION_BED_TO_BIGBED(sorted.sorted)
    validated = VALIDATE_SPLICE_BIGBED(bigbed.bigbed)
    record_input = validated.validated.map { meta, track, source -> tuple(meta + [current_track: 'splice_junction'], track, source) }
    checksummed = CHECKSUM_TRACK(record_input)
    records = WRITE_TRACK_RECORD(checksummed.checksummed)

    emit:
    results = records.results
    versions = records.versions
}
