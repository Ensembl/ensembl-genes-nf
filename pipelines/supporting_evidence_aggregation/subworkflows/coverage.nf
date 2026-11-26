include { VERIFY_SOURCE_CHECKSUM } from '../../track_file_aggregation/modules/verify_source_checksum.nf'
include { VALIDATE_BIGWIG as VALIDATE_SUPPORTING_EVIDENCE_BIGWIG } from '../../../modules/track/validate_bigwig.nf'
include { BUILD_COVERAGE_MATRIX } from '../modules/build_coverage_matrix.nf'

def safe_id(String value) {
    value.replaceAll(/[^A-Za-z0-9_.-]/, '_')
}

workflow BUILD_COVERAGE_MATRICES {
    take:
    plan_tsv

    main:
    source_rows = plan_tsv.splitCsv(header: true, sep: '\t', quote: '"')
        .filter { row -> row.track_type in ['coverage_forward', 'coverage_reverse'] }
        .map { row ->
            def strand = row.track_type == 'coverage_forward' ? 'forward' : 'reverse'
            def meta = [
                id: row.entity_id,
                safe_id: safe_id(row.entity_id),
                safe_source_id: safe_id(row.entity_id),
                source_entity_id: row.entity_id,
                gca_accession: row.gca_accession,
                assembly_release: row.assembly_release,
                track_type: row.track_type,
                coverage_strand: strand,
                expected_sha256: row.sha256,
                normalization_parameters: java.net.URLDecoder.decode(row.normalization_parameters, 'UTF-8'),
                source_path: row.source_path,
                summary_bin_size: params.coverage_bin_size as int,
            ]
            tuple(meta, file(row.track_path, checkIfExists: true))
        }
    checked = VERIFY_SOURCE_CHECKSUM(source_rows)
    validated = VALIDATE_SUPPORTING_EVIDENCE_BIGWIG(checked.checked.map { meta, track ->
        tuple(meta + [id: meta.source_entity_id], track)
    })
    grouped = validated.validated
        .map { meta, track -> tuple("${meta.gca_accession}\t${meta.coverage_strand}", meta, track) }
        .groupTuple()
        .map { _key, metas, tracks ->
            def pairs = metas.withIndex().collect { meta, index -> [meta, tracks[index]] }
                .sort { left, right -> left[0].source_entity_id <=> right[0].source_entity_id }
            def first = pairs[0][0]
            def ordered_metas = pairs.collect { pair -> pair[0] }
            tuple(
                first + [
                    id: "${first.gca_accession}.${first.coverage_strand}",
                    safe_id: safe_id("${first.gca_accession}.${first.coverage_strand}"),
                ],
                pairs.collect { pair -> pair[1] },
                ordered_metas.collect { source_meta -> source_meta.source_entity_id },
            )
        }
    matrices = BUILD_COVERAGE_MATRIX(grouped)

    emit:
    outputs = matrices.matrix.mix(matrices.npz).mix(matrices.descriptor)
    versions = checked.versions
        .mix(validated.versions)
        .mix(matrices.versions)
}
