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


process get_recent_annotations {
    input:
    val annotations_csv

    output:
    stdout

    script:
    """
    python "$projectDir/tasks.py" get_recent_annotations --query_file="$projectDir/get_recent_annotations.sql" --annotations_csv="$annotations_csv"
    """
}


process process_annotation {
    errorStrategy 'ignore'

    input:
    val annotation_database

    output:
    stdout

    script:
    """
    python "$projectDir/tasks.py" process_annotation --annotation_database="$annotation_database"
    """
}


workflow {
    annotations_csv_path = "$workDir/annotations.csv"
    get_recent_annotations(annotations_csv_path)
    /*get_recent_annotations.out.view()*/

    process_annotation(get_recent_annotations.out.splitText())
    process_annotation.out.view()
}
