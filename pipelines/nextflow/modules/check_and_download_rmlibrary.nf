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

process CHECK_AND_DOWNLOAD_REPEATMODELER {
    tag "$gca:repeatmodeler"
    label 'default'
    publishDir "${params.cacheDir}/${gca}/rm_library", mode: 'move'
    afterScript "sleep $params.files_latency"

    input:
    tuple val(url), val(gca)

    output:
    tuple val(gca), path("${gca}.repeatmodeler.fa")

    script:
    """
    set -e
    if wget --spider $url 2>/dev/null; then
        wget -O ${gca}.repeatmodeler.fa $url
    else
        echo "File not found: $url" >&2
        exit 1
    fi
    """
}
