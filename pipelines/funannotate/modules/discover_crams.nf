process DISCOVER_CRAMS {
    tag "${meta.assembly_accession}:${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/curl:7.80.0@sha256:a7e3e10f5b162ac7fdf7fd582056037883e8019611d75a138a7808412161e225'
    containerOptions '--bind /etc:/etc'

    input:
    tuple val(meta), val(cram_dir)

    output:
    tuple val(meta), path('cram_urls.tsv'), emit: urls

    script:
    """
    curl -fL --retry 10 --retry-connrefused --retry-delay 30 --connect-timeout 60 -- '${cram_dir}' \\
        | sed -n 's/.*href="\\([^"]*\\.Aligned\\.SortedByCoordinate\\.out\\.cram\\)".*/\\1/p' \\
        | sort > cram_names.txt
    [ -s cram_names.txt ]
    awk -v dir='${cram_dir}' '{ printf "%d\\t%s%s\\n", NR, dir, \$0 }' cram_names.txt > cram_urls.tsv
    """

    stub:
    """
    printf '1\\tstub.cram\\n' > cram_urls.tsv
    """
}
