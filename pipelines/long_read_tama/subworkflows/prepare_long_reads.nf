nextflow.enable.dsl = 2

include { VALIDATE_APPROVED_LONG_READ_MANIFEST } from '../modules/validate_approved_manifest.nf'
include { FASTQ_DL } from '../modules/fastq_dl.nf'
include { ACQUIRE_LONG_READ_BAM } from '../modules/acquire_bam.nf'
include { RUN_PBCCS } from '../modules/run_pbccs.nf'
include { BAM_TO_FASTQ } from '../modules/bam_to_fastq.nf'
include { VALIDATE_FASTQ } from './validate_fastq.nf'

workflow PREPARE_LONG_READS {
    take:
    approved_source
    validator
    cache_dir

    main:
    // approved_source: path approved manifest
    // reads: tuple val(meta), path(reads.fastq.gz)
    VALIDATE_APPROVED_LONG_READ_MANIFEST(approved_source, validator)
    approved_rows = VALIDATE_APPROVED_LONG_READ_MANIFEST.out.manifest.splitCsv(header: true, sep: '\t')
    approved_rows.branch { row ->
        ont: row.classification == 'ONT_FASTQ'
        ccs_fastq: row.classification == 'PACBIO_CCS_FASTQ'
        processed_fastq: row.classification == 'PACBIO_PROCESSED_FASTQ'
        ccs_bam: row.classification == 'PACBIO_CCS_BAM'
        subread_bam: row.classification == 'PACBIO_SUBREAD_BAM'
    }.set { routes }

    direct_fastq = routes.ont.mix(routes.ccs_fastq).mix(routes.processed_fastq).map { row ->
        tuple([id: row.run_accession, tissue: row.tissue, description: row.description,
               classification: row.classification, proposed_action: row.proposed_action,
               minimap2_preset: row.minimap2_preset, expected_header_representation: row.expected_header_representation],
              row.run_accession, row.selected_artifact_md5, row.selected_artifact_basename, row.selected_artifact_uri)
    }
    bam_inputs = routes.ccs_bam.mix(routes.subread_bam).map { row ->
        tuple([id: row.run_accession, tissue: row.tissue, description: row.description,
               classification: row.classification, proposed_action: row.proposed_action,
               minimap2_preset: row.minimap2_preset, expected_header_representation: row.expected_header_representation],
              row.selected_artifact_md5, row.selected_artifact_basename, row.selected_artifact_uri,
              row.required_artifact_uris ?: '', row.required_artifact_md5s ?: '')
    }
    FASTQ_DL(direct_fastq, cache_dir)
    ACQUIRE_LONG_READ_BAM(bam_inputs, cache_dir)
    ACQUIRE_LONG_READ_BAM.out.artifacts.branch { meta, _files ->
        subreads: meta.classification == 'PACBIO_SUBREAD_BAM'
        ccs: meta.classification == 'PACBIO_CCS_BAM'
    }.set { bam_routes }
    RUN_PBCCS(bam_routes.subreads)
    BAM_TO_FASTQ(bam_routes.ccs.mix(RUN_PBCCS.out.bam))
    canonical_fastq = FASTQ_DL.out.fastq.mix(BAM_TO_FASTQ.out.reads)
    VALIDATE_FASTQ(canonical_fastq, validator)

    emit:
    reads = VALIDATE_FASTQ.out.reads
    reports = VALIDATE_FASTQ.out.report
    molecule_audit = VALIDATE_FASTQ.out.molecule_audit
    versions = VALIDATE_APPROVED_LONG_READ_MANIFEST.out.versions
        .mix(FASTQ_DL.out.versions)
        .mix(ACQUIRE_LONG_READ_BAM.out.versions)
        .mix(RUN_PBCCS.out.versions)
        .mix(BAM_TO_FASTQ.out.versions)
        .mix(VALIDATE_FASTQ.out.versions)
}
