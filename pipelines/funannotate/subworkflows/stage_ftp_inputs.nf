include { STAGE_LOCAL_FILES } from '../modules/stage_local_files.nf'
include { DISCOVER_CRAMS } from '../modules/discover_crams.nf'
include { DOWNLOAD_CRAM } from '../modules/download_cram.nf'
include { ASSEMBLE_STAGED_FILES } from '../modules/assemble_staged_files.nf'

workflow STAGE_FTP_INPUTS {
    take:
    records
    stage_dir

    main:
    STAGE_LOCAL_FILES(records)

    cram_dirs = records
        .filter { meta, _local_files -> meta.files.cram_dir }
        .map { meta, _local_files -> tuple(meta, meta.files.cram_dir) }
    DISCOVER_CRAMS(cram_dirs)

    discovered_crams = DISCOVER_CRAMS.out.urls.flatMap { meta, urls_file ->
        urls_file.text
            .readLines()
            .findAll { line -> line.trim() }
            .collect { line ->
                def fields = line.split('\t', 2)
                tuple(meta, fields[0] as Integer, fields[1], '')
            }
    }
    explicit_crams = records.flatMap { meta, _local_files ->
        meta.files.crams.withIndex().collect { spec, i ->
            tuple(meta, i + 1, spec.url, spec.checksum ?: '')
        }
    }
    cram_specs = explicit_crams.mix(discovered_crams)
    DOWNLOAD_CRAM(cram_specs, stage_dir)

    local_staged = STAGE_LOCAL_FILES.out.staged.map { meta, local_stage ->
        tuple(meta.safe_id, meta, local_stage)
    }
    downloaded_crams = DOWNLOAD_CRAM.out.crams
        .groupTuple()
        .map { meta, cram_files -> tuple(meta.safe_id, cram_files) }
    assemble_inputs = local_staged
        .join(downloaded_crams, by: 0)
        .map { _safe_id, meta, local_stage, cram_files ->
            tuple(meta, local_stage, cram_files)
        }
    ASSEMBLE_STAGED_FILES(assemble_inputs, stage_dir)

    emit:
    staged = ASSEMBLE_STAGED_FILES.out.staged
    checksums = ASSEMBLE_STAGED_FILES.out.checksums
}
