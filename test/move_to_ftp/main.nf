#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MOVE_TO_FTP } from '../../modules/move_to_ftp.nf'

workflow {
    // Create a test file
    test_file = Channel.of([
        [id: 'test_transfer'],
        file("${workflow.launchDir}/test_input.txt")
    ])

    // Create the test input file
    CREATE_TEST_FILE()

    // Move to FTP
    MOVE_TO_FTP(
        CREATE_TEST_FILE.out.test_file,
        params.ftp_destination
    )
}

process CREATE_TEST_FILE {
    output:
    tuple val(meta), path("test_file.txt"), emit: test_file

    script:
    meta = [id: 'test_transfer']
    """
    echo "This is a test file for FTP transfer" > test_file.txt
    echo "Generated at: \$(date)" >> test_file.txt
    echo "Hostname: \$(hostname)" >> test_file.txt
    """
}
