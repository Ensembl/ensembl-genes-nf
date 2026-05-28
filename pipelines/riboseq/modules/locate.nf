process LOCATE {
    tag "${meta.id}"
    label 'process_light'

    input:
    tuple val(meta), val(run)

    output:
    tuple val(meta), path("${run}.collapsed.fa.gz"), emit: collapsed_reads, optional: true
    tuple val(meta), val(run), path("${run}_needs_processing"), emit: needs_processing, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    collapsed_file="${params.collapsed_read_path}/${run}_rpfs.collapsed.fa"
    collapsed_file_1="${params.collapsed_read_path}/${run}_1_rpfs.collapsed.fa.gz"

    if [ -f "\$collapsed_file" ]; then
        gzip -c "\$collapsed_file" > "${run}.collapsed.fa.gz"
        echo "Found collapsed file for $run"
    elif [ -f "\$collapsed_file_1" ]; then
        ln -s "\$collapsed_file_1" "${run}.collapsed.fa.gz"
        echo "Found collapsed _1 file for $run"
    else
        echo "Collapsed file not found for $run. Needs processing."
        touch "${run}_needs_processing"
    fi
    """

    stub:
    """
    touch "${run}.collapsed.fa.gz"
    touch "${run}_needs_processing"
    """
}
