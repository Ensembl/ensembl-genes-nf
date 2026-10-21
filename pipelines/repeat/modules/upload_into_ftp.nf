#!/usr/bin/env nextflow

process UPLOAD_INTO_FTP {

    tag "${meta.gca}:upload_into_ftp"
    label 'ensembl_ftp'

    input:
    tuple val(meta), path(files, arity: '1..*')
    val ftp_base

    output:
    tuple val(meta), path(files), emit: library_out
    path "versions.yml", emit: versions_file

    script:
    def destination = "${params.ftb_local}/${meta.species_name}"

    """
    sudo -u genebuild mkdir -p "${destination}"

    for file in ${files.collect { f -> "\"${f}\"" }.join(' ')}; do
        if [[ -f "\$file" ]]; then
            sudo -u genebuild cp "\$file" "${destination}/"
        fi
    done

    cat > versions.yml <<EOF
"${task.process}":
    upload_repeat_library_into_ftp: "1.0"
EOF
    """
}