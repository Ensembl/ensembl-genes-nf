process DOWNLOAD_REPEATMASKER {
    label 'process_medium'

    tag "${meta.id}"
    publishDir { "${params.outdir}/refseq/${meta.id}/repeatmasker" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(refseq_loaded)

    output:
    tuple val(meta), path('repeatmasker.out.gz'), path('repeatmasker_status.txt'), emit: repeatmasker
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    set +e
    gff-loader \\
        --log-file gff-loader.log \\
        refseq \\
        download-repeatmasker \\
        --base-dir repeatmasker_data \\
        --assembly-acc ${meta.id} \\
        ${args} \\
        > repeatmasker_download.out 2> repeatmasker_download.err
    download_status=\$?
    set -e

    if [[ \$download_status -ne 0 ]]; then
        if grep -Eq '404 Client Error|404 Not Found' repeatmasker_download.err; then
            echo "missing" > repeatmasker_status.txt
            : > repeatmasker.out.gz
            echo "RepeatMasker file not found for ${meta.id}; skipping RepeatMasker loading" >&2
        else
            cat repeatmasker_download.err >&2
            exit \$download_status
        fi
    else
        repeatmasker_file=\$(find repeatmasker_data -type f -name '*_rm.out.gz' -print -quit)
        if [[ -z "\$repeatmasker_file" ]]; then
            echo "RepeatMasker download completed but produced no file for ${meta.id}" >&2
            exit 1
        fi
        cp "\$repeatmasker_file" repeatmasker.out.gz
        echo "available" > repeatmasker_status.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print \$2}')
        ensembl-genes: \$(python -c 'from importlib.metadata import version; print(version("ensembl-genes"))')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p repeatmasker_data/${meta.id}
    touch repeatmasker.out.gz
    echo "available" > repeatmasker_status.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
        ensembl-genes: stub
    END_VERSIONS
    """
}
