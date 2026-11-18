include { STAGE_FTP_INPUTS } from './stage_ftp_inputs.nf'
include { INDEX_AND_VALIDATE_BAMS } from '../modules/index_and_validate_bams.nf'
include { CONVERT_CRAMS } from '../modules/convert_crams.nf'
include { MERGE_BAMS } from '../modules/merge_bams.nf'
include { CONVERT_PROTEIN_ALIGNMENTS } from '../modules/convert_protein_alignments.nf'
include { TRANSFORM_SCALLOP } from '../modules/transform_scallop.nf'
include { RUN_FUNANNOTATE } from '../modules/run_funannotate.nf'

workflow RUN_FUNANNOTATE_PIPELINE {
    take:
    records
    stage_dir
    converter
    cpus

    main:
    // records: channel of [ val(meta) ], with all evidence specifications in meta.files
    STAGE_FTP_INPUTS(records, stage_dir)
    INDEX_AND_VALIDATE_BAMS(STAGE_FTP_INPUTS.out.staged)
    CONVERT_CRAMS(INDEX_AND_VALIDATE_BAMS.out.indexed)
    MERGE_BAMS(CONVERT_CRAMS.out.converted)
    CONVERT_PROTEIN_ALIGNMENTS(
        MERGE_BAMS.out.merged,
        converter
    )
    TRANSFORM_SCALLOP(CONVERT_PROTEIN_ALIGNMENTS.out.converted)
    RUN_FUNANNOTATE(
        TRANSFORM_SCALLOP.out.transformed,
        cpus
    )

    emit:
    results = RUN_FUNANNOTATE.out.results
    gff3 = RUN_FUNANNOTATE.out.gff3
    stdout = RUN_FUNANNOTATE.out.stdout
    stderr = RUN_FUNANNOTATE.out.stderr
    versions = RUN_FUNANNOTATE.out.versions
    checksums = STAGE_FTP_INPUTS.out.checksums
}
