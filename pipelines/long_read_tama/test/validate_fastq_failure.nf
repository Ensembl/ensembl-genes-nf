nextflow.enable.dsl = 2

include { VALIDATE_FASTQ } from '../modules/validate_fastq.nf'

workflow {
    invalid = channel.of(tuple(
        [id: 'invalid-fastq', classification: 'ONT_FASTQ', expected_header_representation: 'ONT'],
        file("${baseDir}/invalid.fastq.gz")
    ))
    VALIDATE_FASTQ(invalid, file("${baseDir}/../bin/read_input_classification.py"))
}
