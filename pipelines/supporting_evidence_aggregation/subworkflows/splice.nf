include { VERIFY_SOURCE_CHECKSUM } from '../../track_file_aggregation/modules/verify_source_checksum.nf'
include { VALIDATE_SPLICE_BIGBED as VALIDATE_SUPPORTING_EVIDENCE_SPLICE } from '../../../modules/track/validate_splice_bigbed.nf'
include { BIGBED_TO_BED } from '../../track_file_aggregation/modules/bigbed_to_bed.nf'
include { BUILD_SPLICE_MATRIX } from '../modules/build_splice_matrix.nf'

def safe_id(String value) {
    value.replaceAll(/[^A-Za-z0-9_.-]/, '_')
}

workflow BUILD_SPLICE_MATRIX_WORKFLOW {
    take:
    plan_tsv
    builder

    main:
    source_rows = plan_tsv.splitCsv(header: true, sep: '\t', quote: '"')
        .filter { row -> row.track_type == 'splice_junction' }
        .map { row ->
            def meta = [
                id: row.entity_id,
                safe_id: safe_id(row.entity_id),
                safe_source_id: safe_id(row.entity_id),
                source_entity_id: row.entity_id,
                gca_accession: row.gca_accession,
                assembly_release: row.assembly_release,
                track_type: row.track_type,
                expected_sha256: row.sha256,
                source_path: row.source_path,
            ]
            tuple(meta, file(row.track_path, checkIfExists: true))
        }
    checked = VERIFY_SOURCE_CHECKSUM(source_rows)
    validated = VALIDATE_SUPPORTING_EVIDENCE_SPLICE(checked.checked)
    beds = BIGBED_TO_BED(validated.validated.map { meta, bigbed ->
        tuple(meta + [aggregate_id: meta.gca_accession], bigbed)
    })
    grouped = beds.bed
        .map { meta, bed -> tuple(meta.gca_accession, meta, bed) }
        .groupTuple()
        .map { _key, metas, input_beds ->
            def pairs = metas.withIndex().collect { meta, index -> [meta, input_beds[index]] }
                .sort { left, right -> left[0].source_entity_id <=> right[0].source_entity_id }
            def first = pairs[0][0]
            tuple(
                first + [
                    id: safe_id(first.gca_accession),
                    safe_id: safe_id(first.gca_accession),
                ],
                pairs.collect { pair -> pair[1] },
                pairs.collect { pair -> pair[0].source_entity_id },
            )
        }
    matrices = BUILD_SPLICE_MATRIX(grouped, builder)

    emit:
    outputs = matrices.matrix.mix(matrices.metadata).mix(matrices.descriptor)
    versions = checked.versions
        .mix(validated.versions)
        .mix(beds.versions)
        .mix(matrices.versions)
}
