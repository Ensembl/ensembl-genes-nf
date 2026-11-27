process DOWNLOAD_CRAM {
    tag "${meta.assembly_accession}:${meta.id}:cram_${cram_index}"
    label 'process_medium'
    container 'quay.io/biocontainers/curl:7.80.0@sha256:a7e3e10f5b162ac7fdf7fd582056037883e8019611d75a138a7808412161e225'
    containerOptions '--bind /etc:/etc'

    input:
    tuple val(meta), val(cram_index), val(cram_url), val(cram_checksum)
    val stage_dir

    output:
    tuple val(meta), path('cram_*.cram'), emit: crams

    script:
    def cram_name = "cram_${String.format('%02d', cram_index as Integer)}.cram"
    def persistent = "${stage_dir}/${meta.safe_id}/${cram_name}"
    def partial = "${persistent}.partial"
    def checksum = cram_checksum ?: ''
    """
    mkdir -p '${stage_dir}/${meta.safe_id}'
    if [ -s '${persistent}' ]; then
        cp -p '${persistent}' '${cram_name}'
    else
        curl -fL --retry 10 --retry-connrefused --retry-delay 30 --connect-timeout 60 -C - -o '${partial}' -- '${cram_url}'
        mv '${partial}' '${persistent}'
        cp -p '${persistent}' '${cram_name}'
    fi
    if [ -n '${checksum}' ]; then
        if ! echo '${checksum}  ${cram_name}' | sha256sum -c -; then
            rm -f '${persistent}' '${cram_name}'
            exit 1
        fi
    fi
    """

    stub:
    """
    touch cram_${String.format('%02d', cram_index as Integer)}.cram
    """
}
