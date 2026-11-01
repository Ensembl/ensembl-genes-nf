include { VERIFY_SOURCE_CHECKSUM } from '../modules/verify_source_checksum.nf'
include { MERGE_BIGWIGS } from '../modules/merge_bigwigs.nf'
include { SCALE_BEDGRAPH } from '../modules/scale_bedgraph.nf'
include { BEDGRAPH_TO_BIGWIG } from '../modules/bedgraph_to_bigwig.nf'
include { VALIDATE_BIGWIG } from '../../../modules/track/validate_bigwig.nf'
include { CHECKSUM_TRACK } from '../../../modules/track/checksum_track.nf'
include { WRITE_AGGREGATE_RECORD } from '../modules/write_aggregate_record.nf'

def safe_entity_id(String value) {
    value.replaceAll(/[^A-Za-z0-9_.-]/, '_')
}

workflow AGGREGATE_COVERAGE_TRACKS {
    take:
    plan_tsv // path: one expanded source-track row per line
    chrom_sizes // singleton path: one assembly reference
    record_script // singleton path: aggregate record helper

    main:
    source_rows = plan_tsv.splitCsv(header: true, sep: '\t', quote: '"')
        .filter { row -> row.track_type in ['coverage_forward', 'coverage_reverse'] }
        .map { row ->
            def strand = row.track_type == 'coverage_forward' ? 'forward' : 'reverse'
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
                coverage_strand: strand,
                expected_sha256: row.expected_sha256,
                normalization_parameters: java.net.URLDecoder.decode(row.normalization_parameters, 'UTF-8'),
                normalization_signature: java.net.URLDecoder.decode(row.normalization_signature, 'UTF-8'),
                aggregation_method: row.aggregation_method,
                metadata_json: java.net.URLDecoder.decode(row.metadata_json, 'UTF-8'),
                source_path: row.source_path,
            ]
            tuple(meta, file(row.track_path, checkIfExists: true))
        }
    checked = VERIFY_SOURCE_CHECKSUM(source_rows)
    validated = VALIDATE_BIGWIG(checked.checked.map { meta, track -> tuple(meta, track) })
    grouped = validated.validated
        .map { meta, track -> tuple("${meta.gca_accession}\t${meta.aggregate_id}\t${meta.coverage_strand}", meta, track) }
        .groupTuple()
        .map { _key, metas, tracks ->
            def first = metas[0]
            def source_ids = metas.collect { source_meta -> source_meta.source_entity_id }.unique().sort()
            def source_paths = metas.collect { source_meta -> source_meta.source_path }.unique().sort()
            def normalization = [
                aggregation_method: first.aggregation_method,
                source_normalization_parameters: first.normalization_parameters,
                source_count: tracks.size(),
                strand: first.coverage_strand,
            ]
            try {
                def source_parameters = new groovy.json.JsonSlurper().parseText(first.normalization_parameters)
                if (source_parameters instanceof Map) {
                    ['bin_size', 'normalization'].each { key ->
                        if (source_parameters.containsKey(key)) {
                            normalization[key] = source_parameters[key]
                        }
                    }
                }
            } catch (Exception _ignored) {
                // Legacy manifests may contain opaque normalization metadata.
                // The validator records that compatibility mode explicitly.
            }
            def group_meta = first + [
                id: first.aggregate_id,
                current_track: "coverage_${first.coverage_strand}",
                output_name: "${safe_entity_id(first.aggregate_id)}.${first.coverage_strand}.bw",
                source_count: tracks.size(),
                source_entity_ids: source_ids.join(','),
                source_path: source_paths.join(','),
                normalization_parameters: groovy.json.JsonOutput.toJson(normalization),
                weight_divisor: first.aggregation_method == 'equal_sample' ? tracks.size() : 1,
            ]
            tuple(group_meta, tracks)
        }
    merged = MERGE_BIGWIGS(grouped)
    raw_bedgraphs = merged.bedgraph.filter { meta, _bedgraph -> meta.aggregation_method == 'raw' }
    weighted_bedgraphs = merged.bedgraph.filter { meta, _bedgraph -> meta.aggregation_method == 'equal_sample' }
    scaled = SCALE_BEDGRAPH(weighted_bedgraphs)
    bedgraphs = raw_bedgraphs.mix(scaled.scaled)
    bigwig_input = bedgraphs.map { meta, bedgraph -> tuple(meta, bedgraph, chrom_sizes) }
    bigwigs = BEDGRAPH_TO_BIGWIG(bigwig_input)
    output_meta = bigwigs.bigwig.map { meta, track -> tuple(meta + [current_track: "coverage_${meta.coverage_strand}"], track) }
    checksummed = CHECKSUM_TRACK(output_meta)
    records = WRITE_AGGREGATE_RECORD(checksummed.checksummed, record_script)

    emit:
    results = records.results
    versions = checked.versions
        .mix(merged.versions)
        .mix(scaled.versions)
        .mix(bigwigs.versions)
        .mix(checksummed.versions)
        .mix(records.versions)
}
