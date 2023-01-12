#!/usr/bin/env nextflow


/*
See the NOTICE file distributed with this work for additional information
regarding copyright ownership.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/


nextflow.enable.dsl=2


log.info """\
Genebuild annotation statistics Nextflow pipeline

Usage:
    nextflow run annotation_stats.nf --message MESSAGE
"""


// default pipeline parameters
params.message = "Hello World!"


// default workflow
workflow {
    // input data channel
    input_channel = Channel.value(params.message)

    print_message(input_channel)

    // print process output to the terminal
    print_message.out.view()
}


// https://www.nextflow.io/docs/latest/process.html
process print_message {
    input:
    val message

    output:
    stdout

    script:
    """
    echo "${params.message}"
    """
}
