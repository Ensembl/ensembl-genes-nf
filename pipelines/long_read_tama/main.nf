#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VALIDATE_LONG_READ_MANIFEST } from './modules/validate_manifest.nf'
include { DISCOVER_LONG_READ_DATA } from './modules/discover_long_read.nf'
include { INSPECT_LONG_READ_MANIFEST } from './modules/inspect_manifest.nf'
include { RESOLVE_LONG_READ_METADATA } from './modules/resolve_metadata.nf'
include { VALIDATE_APPROVED_LONG_READ_MANIFEST } from './modules/validate_approved_manifest.nf'
include { FASTQ_DL } from './modules/fastq_dl.nf'
include { ACQUIRE_LONG_READ_BAM } from './modules/acquire_bam.nf'
include { CANONICALISE_PACBIO_READS } from './modules/canonicalise_pacbio.nf'
include { PREPARE_LONG_READ_INPUT } from './modules/prepare_fastq.nf'
include { BUILD_MINIMAP2_INDEX } from './modules/build_minimap_index.nf'
include { MINIMAP2_TO_SORTED_BAM } from './modules/minimap2_sorted_bam.nf'
include { SPLIT_BAM_BY_CONTIG } from './modules/split_contigs.nf'
include { TAMA_COLLAPSE } from './modules/tama_collapse.nf'
include { TAMA_MERGE } from './modules/tama_merge.nf'
include { TAMA_MERGE as TAMA_MERGE_ACCESSION } from './modules/tama_merge.nf'
include { VALIDATE_LONG_READ_MODELS } from './modules/validate_models.nf'

def make_read_meta(row, cohort) {
    [id: row.run_accession, tissue: row.tissue, description: row.description,
     classification: row.classification, proposed_action: row.proposed_action,
     minimap2_preset: row.minimap2_preset, expected_header_representation: row.expected_header_representation,
     cohort: cohort]
}

workflow {
    if (!params.manifest && !params.approved_manifest && !params.taxon_id) error 'Provide --taxon_id for discovery, --manifest for inventory, or --approved_manifest for production processing'
    if (params.approved_manifest && !params.reference_fasta) error 'Provide --reference_fasta with the reference FASTA for production processing'
    if (!(params.shard_mode in ['none', 'contig'])) error "params.shard_mode must be 'none' or 'contig'"
    if (!(params.secondary_mode in ['no', 'yes'])) error "params.secondary_mode must be 'no' or 'yes'"
    if (!params.fastq_cache_dir && params.approved_manifest) error 'Provide --fastq_cache_dir for production processing'

    validator = file("${baseDir}/bin/read_input_classification.py", checkIfExists: true)
    if (params.taxon_id && !params.manifest && !params.approved_manifest) {
        DISCOVER_LONG_READ_DATA(params.taxon_id, file("${baseDir}/bin/discover_long_read_data.py", checkIfExists: true), params.discovery_cache_dir ?: "${params.outdir}/discovery_cache")
        return
    }
    if (!params.approved_manifest) {
        VALIDATE_LONG_READ_MANIFEST(file(params.manifest, checkIfExists: true), file("${baseDir}/bin/normalise_manifest.py", checkIfExists: true))
        if (params.metadata_json) {
            inventory_metadata = file(params.metadata_json, checkIfExists: true)
        } else {
            metadata_cache = params.metadata_cache_dir ?: "${params.outdir}/metadata_cache"
            RESOLVE_LONG_READ_METADATA(
                VALIDATE_LONG_READ_MANIFEST.out.manifest,
                file("${baseDir}/bin/resolve_metadata.py", checkIfExists: true),
                metadata_cache
            )
            inventory_metadata = RESOLVE_LONG_READ_METADATA.out.metadata
        }
        INSPECT_LONG_READ_MANIFEST(VALIDATE_LONG_READ_MANIFEST.out.manifest, inventory_metadata, validator)
        return
    }

    reference_fasta = file(params.reference_fasta, checkIfExists: true)
    approved_source = file(params.approved_manifest, checkIfExists: true)
    VALIDATE_APPROVED_LONG_READ_MANIFEST(approved_source, validator)
    approved_rows = VALIDATE_APPROVED_LONG_READ_MANIFEST.out.manifest.splitCsv(header: true, sep: '\t')
    approved_rows.branch { row ->
        ont: row.classification == 'ONT_FASTQ'
        ccs_fastq: row.classification == 'PACBIO_CCS_FASTQ'
        processed_fastq: row.classification == 'PACBIO_PROCESSED_FASTQ'
        ccs_bam: row.classification == 'PACBIO_CCS_BAM'
        subread_bam: row.classification == 'PACBIO_SUBREAD_BAM'
    }.set { routes }
    direct_fastq = routes.ont.mix(routes.ccs_fastq).mix(routes.processed_fastq).map { row -> tuple(make_read_meta(row, params.cohort_id), row.run_accession, row.selected_artifact_md5, row.selected_artifact_basename, row.selected_artifact_uri) }
    bam_inputs = routes.ccs_bam.mix(routes.subread_bam).map { row -> tuple(make_read_meta(row, params.cohort_id), row.selected_artifact_md5, row.selected_artifact_basename, row.selected_artifact_uri, row.required_artifact_uris ?: '', row.required_artifact_md5s ?: '') }
    cache = file(params.fastq_cache_dir).toAbsolutePath().toString()
    FASTQ_DL(direct_fastq, cache)
    ACQUIRE_LONG_READ_BAM(bam_inputs, cache)
    CANONICALISE_PACBIO_READS(ACQUIRE_LONG_READ_BAM.out.artifacts.map { meta, files -> tuple(meta, files) }, validator)
    PREPARE_LONG_READ_INPUT(FASTQ_DL.out.fastq.mix(CANONICALISE_PACBIO_READS.out.reads), validator)

    if (params.minimap_index) minimap_index = file(params.minimap_index, checkIfExists: true)
    else { BUILD_MINIMAP2_INDEX(reference_fasta); minimap_index = BUILD_MINIMAP2_INDEX.out.index }
    MINIMAP2_TO_SORTED_BAM(PREPARE_LONG_READ_INPUT.out.reads, minimap_index, reference_fasta)
    if (params.shard_mode == 'contig') {
        SPLIT_BAM_BY_CONTIG(MINIMAP2_TO_SORTED_BAM.out.bam)
        contig_bams = SPLIT_BAM_BY_CONTIG.out.shards.flatMap { meta, shard_dir -> shard_dir.listFiles().findAll { bam -> bam.name.endsWith('.bam') }.collect { bam -> tuple(meta, bam.baseName, bam) } }
        TAMA_COLLAPSE(contig_bams, reference_fasta)
        accession_beds = TAMA_COLLAPSE.out.bed.map { meta, _shard, bed -> tuple(meta.id, bed) }.groupTuple().map { accession, beds -> tuple(accession, beds) }
        TAMA_MERGE_ACCESSION(accession_beds)
        final_beds = TAMA_MERGE_ACCESSION.out.bed.collect()
    } else {
        whole_bams = MINIMAP2_TO_SORTED_BAM.out.bam.map { meta, bam, _bai -> tuple(meta, 'whole', bam) }
        TAMA_COLLAPSE(whole_bams, reference_fasta)
        final_beds = TAMA_COLLAPSE.out.bed.map { _meta, _shard, bed -> bed }.collect()
    }
    TAMA_MERGE(final_beds.map { beds -> tuple(params.cohort_id, beds) })
    VALIDATE_LONG_READ_MODELS(TAMA_MERGE.out.bed)
}
