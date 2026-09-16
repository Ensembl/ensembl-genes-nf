nextflow.enable.dsl = 2

include { STATS_FASTQ } from '../modules/stats_fastq.nf'
include { PROBE_FASTQ } from '../modules/probe_fastq.nf'
include { VALIDATE_FASTQ as VALIDATE_FASTQ_RECORDS } from '../modules/validate_fastq_records.nf'
include { AUDIT_FASTQ } from '../modules/audit_fastq.nf'
include { STAGE_VALIDATED_FASTQ } from '../modules/stage_validated_fastq.nf'

workflow VALIDATE_FASTQ {
    take:
    reads
    validator

    main:
    STATS_FASTQ(reads, validator)
    PROBE_FASTQ(reads, validator)
    VALIDATE_FASTQ_RECORDS(reads, validator)

    audit_input = PROBE_FASTQ.out.probe.join(VALIDATE_FASTQ_RECORDS.out.validation, by: 0)
    AUDIT_FASTQ(audit_input)

    stage_input = reads.join(VALIDATE_FASTQ_RECORDS.out.validation, by: 0)
    STAGE_VALIDATED_FASTQ(stage_input)

    emit:
    reads = STAGE_VALIDATED_FASTQ.out.reads
    report = STATS_FASTQ.out.report
    molecule_audit = AUDIT_FASTQ.out.molecule_audit
    versions = STATS_FASTQ.out.versions
        .mix(PROBE_FASTQ.out.versions)
        .mix(VALIDATE_FASTQ_RECORDS.out.versions)
        .mix(AUDIT_FASTQ.out.versions)
        .mix(STAGE_VALIDATED_FASTQ.out.versions)
}
