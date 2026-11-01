include { VERIFY_SOURCE_CHECKSUM } from '../modules/verify_source_checksum.nf'
include { BIGBED_TO_BED } from '../modules/bigbed_to_bed.nf'
include { AGGREGATE_JUNCTIONS } from '../modules/aggregate_junctions.nf'
include { JUNCTION_BED_TO_BIGBED } from '../../../modules/track/junction_bed_to_bigbed.nf'
include { VALIDATE_SPLICE_BIGBED; VALIDATE_SPLICE_BIGBED as VALIDATE_AGGREGATED_SPLICE_BIGBED } from '../../../modules/track/validate_splice_bigbed.nf'
include { CHECKSUM_TRACK } from '../../../modules/track/checksum_track.nf'
include { WRITE_AGGREGATE_RECORD } from '../modules/write_aggregate_record.nf'

def safe_entity_id(String value) {
    value.replaceAll(/[^A-Za-z0-9_.-]/, '_')
}

workflow AGGREGATE_SPLICE_TRACKS {
    take:
    plan_tsv // path: one expanded source-track row per line
    chrom_sizes // singleton path: one assembly reference
    autosql // singleton path: shared splice-junction AutoSql
    aggregator // singleton path: junction reduction script
    record_script // singleton path: aggregate record helper

    main:
    source_rows = plan_tsv.splitCsv(header: true, sep: '\t', quote: '"')
        .filter { row -> row.track_type == 'splice_junction' }
        .map { row ->
            def meta = [
                id: row.source_entity_id,
                safe_source_id: safe_entity_id(row.source_entity_id),
                source_entity_id: row.source_entity_id,
                source_entity_type: row.source_entity_type,
                gca_accession: row.gca_accession,
                assembly_release: row.assembly_release,
                aggregate_id: row.aggregate_id,
                aggregate_type: row.aggregate_type,
                aggregate_label: row.aggregate_label,
                safe_id: safe_entity_id(row.aggregate_id),
                track_type: row.track_type,
                expected_sha256: row.expected_sha256,
                normalization_parameters: 'not_applicable',
                normalization_signature: 'not_applicable',
                aggregation_method: row.aggregation_method,
                metadata_json: java.net.URLDecoder.decode(row.metadata_json, 'UTF-8'),
                source_path: row.source_path,
            ]
            tuple(meta, file(row.track_path, checkIfExists: true))
        }
    checked = VERIFY_SOURCE_CHECKSUM(source_rows)
    validated = VALIDATE_SPLICE_BIGBED(checked.checked.map { meta, track -> tuple(meta, track) })
    beds = BIGBED_TO_BED(validated.validated)
    grouped = beds.bed
        .map { meta, bed -> tuple("${meta.gca_accession}\t${meta.aggregate_id}", meta, bed) }
        .groupTuple()
        .map { _key, metas, input_beds ->
            def first = metas[0]
            def source_ids = metas.collect { source_meta -> source_meta.source_entity_id }.unique().sort()
            def source_paths = metas.collect { source_meta -> source_meta.source_path }.unique().sort()
            def group_meta = first + [
                id: first.aggregate_id,
                current_track: 'splice_junction',
                output_name: "${safe_entity_id(first.aggregate_id)}.bb",
                source_count: source_ids.size(),
                source_entity_ids: source_ids.join(','),
                source_path: source_paths.join(','),
                normalization_parameters: 'not_applicable',
            ]
            tuple(group_meta, input_beds)
        }
    aggregated = AGGREGATE_JUNCTIONS(grouped, aggregator)
    bigbed_input = aggregated.bed.map { meta, bed -> tuple(meta, bed, chrom_sizes) }
    bigbeds = JUNCTION_BED_TO_BIGBED(bigbed_input, autosql)
    output_meta = bigbeds.bigbed.map { meta, track -> tuple(meta + [current_track: 'splice_junction'], track) }
    validated_output = VALIDATE_AGGREGATED_SPLICE_BIGBED(output_meta)
    checksummed = CHECKSUM_TRACK(validated_output.validated)
    records = WRITE_AGGREGATE_RECORD(checksummed.checksummed, record_script)

    emit:
    results = records.results
    versions = checked.versions
        .mix(beds.versions)
        .mix(aggregated.versions)
        .mix(bigbeds.versions)
        .mix(checksummed.versions)
        .mix(records.versions)
}
