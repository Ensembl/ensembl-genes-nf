#!/usr/bin/env nextflow

/*
@category Check_and download RepeatModeler library

@summary
Check if a RepeatModeler library file exists at the given URL and download it if it does.

@description
This process checks if a RepeatModeler library file exists at the given URL.
If the file exists, it downloads the file and saves it with the name "<GCA>.
repeatmodeler.fa". If the file does not exist, it outputs an error message 
and exits with a non-zero status.

@implementation
- Validates the remote URL using `wget --spider`.
- Downloads the library as `<GCA>.repeatmodeler.fa`.
- Stops the workflow if the download target is unavailable.
- Generates `versions.yml` containing the wget version.

@inputs
meta             Genome metadata map
url              RepeatModeler library URL

@outputs
meta          Genome metadata map
fasta_file   RepeatModeler library file
*/
process CHECK_AND_DOWNLOAD_RMLIBRARY {
    tag "$meta.gca"
    label 'default'
    publishDir "${params.outdir}/${meta.gca}/rm_library", mode: 'copy'
    afterScript "sleep $params.files_latency"

    input:
    tuple val(url), val(meta)

    output:
    tuple val(meta), path("*.fa"), emit: repeatmodeler_library_out
    path "versions.yml", emit: versions_file

    script:
    """
    set -e
    if wget --spider $url 2>/dev/null; then
        wget -O ${meta.gca}.repeatmodeler.fa $url
    else
        echo "File not found: $url" >&2
        exit 1
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version 2>&1 | head -n 1 | sed 's/GNU Wget //')
    END_VERSIONS
    """
}
