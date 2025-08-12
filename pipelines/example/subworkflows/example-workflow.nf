include { DUMMY_PROCESS } from './processes/dummy_process.nf'

workflow MY_SUBWORKFLOW {

    take:
    input_channel

    main:
    DUMMY_PROCESS(input_channel)

    emit:
    results = DUMMY_PROCESS.out.results
}