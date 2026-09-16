nextflow.enable.dsl = 2

include { VALIDATE_FASTQ } from '../subworkflows/validate_fastq.nf'

process MAKE_TINY_FASTQ {
    input:
    val meta

    output:
    tuple val(meta), path('tiny.fastq.gz'), emit: reads

    script:
    """
    printf '@123e4567-e89b-12d3-a456-426614174000 runid=abc ch=4\\nACGT\\n+\\n!!!!\\n@123e4567-e89b-12d3-a456-426614174001 runid=abc ch=4\\nTGCA\\n+\\n!!!!\\n' | gzip -c > tiny.fastq.gz
    """
}

workflow {
    meta = [id: 'tiny', classification: 'ONT_FASTQ', expected_header_representation: 'ONT']
    MAKE_TINY_FASTQ(meta)
    VALIDATE_FASTQ(
        MAKE_TINY_FASTQ.out.reads,
        file("${baseDir}/../bin/read_input_classification.py", checkIfExists: true)
    )
}
