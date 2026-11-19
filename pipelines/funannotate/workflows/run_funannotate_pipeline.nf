include { PREPARE_FUNANNOTATE_INPUTS } from '../subworkflows/prepare_funannotate_inputs.nf'
include { RUN_FUNANNOTATE } from '../modules/run_funannotate.nf'

workflow RUN_FUNANNOTATE_PIPELINE {
    take:
    records
    stage_dir
    cpus

    main:
    PREPARE_FUNANNOTATE_INPUTS(records, stage_dir)
    RUN_FUNANNOTATE(
        PREPARE_FUNANNOTATE_INPUTS.out.prepared,
        cpus,
    )

    emit:
    results = RUN_FUNANNOTATE.out.results
    gff3 = RUN_FUNANNOTATE.out.gff3
    stdout = RUN_FUNANNOTATE.out.stdout
    stderr = RUN_FUNANNOTATE.out.stderr
    versions = RUN_FUNANNOTATE.out.versions
    checksums = PREPARE_FUNANNOTATE_INPUTS.out.checksums
}
