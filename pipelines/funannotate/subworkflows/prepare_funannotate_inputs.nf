include { STAGE_FTP_INPUTS } from './stage_ftp_inputs.nf'
include { INDEX_AND_VALIDATE_BAMS } from '../modules/index_and_validate_bams.nf'
include { CONVERT_CRAMS } from '../modules/convert_crams.nf'
include { MERGE_BAMS } from '../modules/merge_bams.nf'
include { CONVERT_PROTEIN_ALIGNMENTS } from '../modules/convert_protein_alignments.nf'
include { TRANSFORM_SCALLOP } from '../modules/transform_scallop.nf'

workflow PREPARE_FUNANNOTATE_INPUTS {
    take:
    records
    stage_dir

    main:
    // records: channel of tuple(meta, local_files), with evidence specifications in meta.files
    STAGE_FTP_INPUTS(records, stage_dir)
    INDEX_AND_VALIDATE_BAMS(STAGE_FTP_INPUTS.out.staged)
    CONVERT_CRAMS(INDEX_AND_VALIDATE_BAMS.out.indexed)
    MERGE_BAMS(CONVERT_CRAMS.out.converted)
    CONVERT_PROTEIN_ALIGNMENTS(MERGE_BAMS.out.merged)
    TRANSFORM_SCALLOP(CONVERT_PROTEIN_ALIGNMENTS.out.converted)

    emit:
    prepared = TRANSFORM_SCALLOP.out.transformed
    checksums = STAGE_FTP_INPUTS.out.checksums
}
