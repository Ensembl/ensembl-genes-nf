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
    nextflow run annotation_stats.nf [...]
"""


process check_stats_files {
    input:
    val annotation_directory
    val production_name

    output:
    stdout

    script:
    """
    python "$projectDir/tasks.py" check_stats_files "$annotation_directory" "$production_name"
    """
}


workflow {
    annotation_directory = Channel.fromPath(params.annotation_directory)
    production_name = Channel.value(params.production_name)
    check_stats_files(annotation_directory, production_name)

    check_stats_files.out.view()
}
