#!/usr/bin/env nextflow

/*
This process fetches the RepeatModeler library file for a given GCA accession. 
If the file is available at the specified URL, it downloads the file 
and saves it with the name "<GCA>.repeatmodeler.fa". 
If the file is not available, it outputs a message indicating that 
the download was skipped and creates an empty file with the same name. 
The fetched library file is saved in the "library" directory under 
the output directory for the given GCA accession.
*/

process CHECK_FTP {
    tag "${meta.gca}:check_${file_type}"
    label 'fetch_file'

    input:
    tuple val(meta), val(url), val(file_type)

    output:
    tuple val(meta), val(url), val(file_type), path("status.txt"),
        emit: status_out

    path "versions.yml", emit: versions_file

    script:
    """
    HTTP_STATUS=\$(curl \\
        --silent \\
        --show-error \\
        --location \\
        --output /dev/null \\
        --write-out '%{http_code}' \\
        --head \\
        --retry 5 \
    --retry-all-errors \
    --retry-delay 10 \
    --connect-timeout 15 \
    --max-time 60 \
        "${url}")

    if [ "\$HTTP_STATUS" = "200" ]; then
        echo "true" > status.txt
        echo "${file_type} exists for ${meta.gca}"

    elif [ "\$HTTP_STATUS" = "404" ]; then
        echo "false" > status.txt
        echo "${file_type} missing for ${meta.gca}"

    else
        echo "Unexpected HTTP status: \$HTTP_STATUS" >&2
        exit 1
    fi

    CURL_VERSION=\$(curl --version | head -n 1 | awk '{print \$2}')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: \$CURL_VERSION
    END_VERSIONS
    """
    }

