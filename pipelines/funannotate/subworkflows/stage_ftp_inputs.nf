include { STAGE_FTP_FILES } from '../modules/stage_ftp_files.nf'

workflow STAGE_FTP_INPUTS {
    take:
    records
    stage_dir

    main:
    STAGE_FTP_FILES(records, stage_dir)

    emit:
    staged = STAGE_FTP_FILES.out.staged
    checksums = STAGE_FTP_FILES.out.checksums
}
